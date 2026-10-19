## TSSs of the + and - strand transcripts at intergenic enhancer DHSs from nascent RNA
## (GRO-seq or 4sU-seq), and fraction of enhancers whose two transcripts switch from
## convergent (TSS+ upstream of TSS-) to divergent upon differentiation.
##
## Usage (from 07-Enhancer/Enhancer-analysis/): Rscript 03_nascent_tss_edges.r <dataset>
##   dataset : GSE142312 (GRO-seq, ESC vs NPC) or GSE115774 (4sU-seq, Day0 vs Day3)
##
## Enhancer DHSs not overlapping a protein-coding / lncRNA gene or its promoter (-2000/+500 bp).
## Replicates of each condition are pooled. Per DHS and strand, the TSS is the 5' edge of the
## most 5' block of signal >= 10% of the maximum in the DHS (>= 20 nt long, with >= 33% of the
## signal of the strongest block); DHSs whose signal crosses the 5' boundary are discarded.
## Switch: TSS+ - TSS- < 0 in the control condition and > 0 after differentiation, among DHSs
## with coverage > min_cov on both strands in both conditions (threshold set per dataset).
##
## Input : nascent RNA bigWigs (mm39), ENHANCER_ANNOTATION_RDS, GENCODE_GTF (config/paths.R)
## Output (in $ENHANCER_DIR/nascent_tss_edges/<dataset>/):
##   <condition>_plus.bw, <condition>_minus.bw  pooled coverage at the enhancer DHSs
##   <condition>_TSS_edges.bedGraph             TSS+ (score 1) and TSS- (score -1) per DHS
##   tss_edges_per_DHS.tsv                      TSSs and coverage per DHS and condition
##   enhancer_tss_shift.tsv                     DHSs passing the coverage filter, TSS distance in
##                                              both conditions, switch
##   enhancer_tss_shift_summary.tsv             number and percentage of switching enhancers

library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(parallel)
library(tidyverse)

source("../../config/paths.R")

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]

#### Dataset settings (coverage thresholds set manually per experiment)

settings <- list(
    ## GRO-seq: GSM4225082_ESC_line1_Exp1.pos.bw.mm39.bw.bw, GSM4225085_ESC_NPC_line1_Exp1.neg...
    GSE142312 = list(
        bw_dir = file.path(DATA_DIR, "GSE142312_nascent/mm39"),
        condition = function(b) str_extract(b %>% str_replace("ESC_NPC", "NPC"), "(NPC_line._Exp.|ESC_line._Exp.)"),
        strand = function(b) ifelse(grepl("pos", b), "+", "-"),
        base_condition = function(cond) sub("_.*", "", cond),
        control = "ESC", treated = "NPC", min_cov = 100),
    ## 4sU-seq: GSM3189531_Day0_R1_PosStr.bigWig.bw
    GSE115774 = list(
        bw_dir = file.path(DATA_DIR, "GSE115774_nascent/bw39"),
        condition = function(b) str_extract(b, "Day[0-9]+_R[1-2]"),
        strand = function(b) ifelse(grepl("PosStr", b), "+", "-"),
        base_condition = function(cond) sub("_R[0-9]+$", "", cond),
        control = "Day0", treated = "Day3", min_cov = 200)
)
set <- settings[[dataset]]
if (is.null(set)) stop("dataset must be one of: ", paste(names(settings), collapse = ", "))

out_dir <- file.path(ENHANCER_DIR, "nascent_tss_edges", dataset)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#### Intergenic enhancer DHSs

enhancer.df <- readRDS(ENHANCER_ANNOTATION_RDS)

dhs.list <- enhancer.df %>%
    .$DHS

dhs_gr <- dhs.list %>%
    data.frame(DHS = .) %>%
    distinct() %>%
    separate(DHS, into = c("seqnames", "start", "end"), sep = "[_-]", remove = FALSE) %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    plyranges::as_granges()

## exclude enhancers overlapping protein-coding / lncRNA genes and their promoters
gtf <- import(GENCODE_GTF, format = "gtf")
gencode.gene <- gtf[gtf$type == "gene" & gtf$gene_type %in% c("protein_coding", "lncRNA")]
gencode.promoters <- promoters(gencode.gene, upstream = 2000, downstream = 500)
exclusion_mask <- GenomicRanges::reduce(c(gencode.gene, gencode.promoters))

dhs_intergenic <- subsetByOverlaps(dhs_gr, exclusion_mask, invert = TRUE)


#### bigWig files: plus and minus track of each replicate

bw_files <- list.files(set$bw_dir, pattern = "\\.bw$", full.names = TRUE)

meta_nascent <- data.frame(file_path = bw_files) %>%
    mutate(basename = basename(file_path)) %>%
    mutate(condition = set$condition(basename)) %>%
    mutate(strand = set$strand(basename)) %>%
    pivot_wider(id_cols = condition, names_from = strand, values_from = file_path) %>%
    dplyr::rename(plus_bw = `+`, minus_bw = `-`) %>%
    drop_na()

## replicates are pooled per condition (e.g. "Day0_R1" -> "Day0"); only the 2 compared conditions
meta_nascent$base_condition <- set$base_condition(meta_nascent$condition)
meta_nascent <- meta_nascent %>% dplyr::filter(base_condition %in% c(set$control, set$treated))
print(meta_nascent)


#### TSS edge detection

## 5' edge of the main signal block of one strand in one DHS
get_robust_edge <- function(pos, score, region_start, region_end, strand = "+",
                            bg_ratio = 0.1, min_width = 20) {
    ## Define a consistent empty return template to keep data.table happy
    null_res <- list(tss = NA_integer_, cov = 0)
    if (length(pos) == 0) return(null_res)
    ## 1. Create a continuous, 0-filled array for the exact DHS width
    region_length <- region_end - region_start + 1
    full_cov <- numeric(region_length)
    ## Map the sparse BigWig positions into the continuous array
    valid_idx <- pos - region_start + 1
    in_bounds <- valid_idx >= 1 & valid_idx <= region_length
    full_cov[valid_idx[in_bounds]] <- score[in_bounds]
    ## Strict edge purge: discard the DHS if the signal crosses its 5' boundary (bleed-in)
    thresh <- max(full_cov) * bg_ratio
    n <- length(full_cov)
    if (strand == "+") {
        if (full_cov[1] > thresh) return(null_res)
    } else if (strand == "-") {
        if (full_cov[n] > thresh) return(null_res)
    }
    ## If there is no signal overall, exit early
    if (sum(full_cov) == 0) return(null_res)
    ## 2. Dynamic threshold: signal must exceed 10% (bg_ratio) of the max peak
    active <- full_cov >= thresh
    ## 3. Find continuous blocks (transcription bubbles)
    runs <- rle(active)
    valid_runs <- which(runs$values == TRUE & runs$lengths >= min_width)
    ## If no valid continuous block is found, return the integral, but no TSS
    if (length(valid_runs) == 0) return(list(tss = NA_integer_, cov = sum(full_cov)))
    ## 4. Find the main transcript body (filtering by block mass)
    run_ends_idx <- cumsum(runs$lengths)
    run_starts_idx <- c(1, run_ends_idx[-length(run_ends_idx)] + 1)
    ## A. Total integral of every valid block
    run_integrals <- sapply(valid_runs, function(i) {
        sum(full_cov[run_starts_idx[i]:run_ends_idx[i]])
    })
    ## B. A block is robust if it has >= 33% of the mass of the strongest block
    run_ratio <- 0.33
    max_integral <- max(run_integrals)
    robust_idx <- which(run_integrals >= (max_integral * run_ratio))
    robust_runs <- valid_runs[robust_idx]
    ## C. Pick the most 5' block among the robust blocks
    if (strand == "+") {
        best_run <- robust_runs[1]
    } else {
        best_run <- robust_runs[length(robust_runs)]
    }
    best_idx_start <- run_starts_idx[best_run]
    best_idx_end <- run_ends_idx[best_run]
    ## 5. Return the leading edge coordinate and the coverage integral
    if (strand == "+") {
        best_tss <- as.integer(region_start + best_idx_start - 1)
    } else {
        best_tss <- as.integer(region_start + best_idx_end - 1)
    }
    return(list(tss = best_tss, cov = sum(full_cov)))
}

## pooled coverage of several bigWigs at the DHSs (minus strand tracks may be negative)
pool_bigwigs <- function(bw_paths, dhs_gr) {
    lapply(bw_paths, import.bw, selection = BigWigSelection(dhs_gr)) %>%
        as("GRangesList") %>%
        unlist() %>%
        plyranges::mutate(score = abs(score)) %>%
        coverage(weight = "score") %>%
        as("GRanges") %>%
        .[.$score > 0]
}

## per-base signal of one strand in each DHS
expand_hits <- function(dhs_gr, gr) {
    hits <- findOverlaps(dhs_gr, gr)
    dt <- data.table(
        DHS_id = queryHits(hits),
        start = start(gr)[subjectHits(hits)],
        end = end(gr)[subjectHits(hits)],
        score = score(gr)[subjectHits(hits)]
    )
    dt_1bp <- dt[, .(pos = start:end, score = score), by = .(DHS_id, block_id = 1:nrow(dt))]
    setkey(dt_1bp, DHS_id)
    dt_1bp
}

calculate_nascent_edges <- function(plus_bw, minus_bw, dhs_gr, condition_name) {
    gr_p <- pool_bigwigs(plus_bw, dhs_gr)
    export.bw(gr_p, file.path(out_dir, paste0(condition_name, "_plus.bw")))
    gr_m <- pool_bigwigs(minus_bw, dhs_gr)
    export.bw(gr_m, file.path(out_dir, paste0(condition_name, "_minus.bw")))

    dt_dhs <- as.data.table(dhs_gr)
    dt_dhs[, DHS_id := .I]
    dt_p_1bp <- expand_hits(dhs_gr, gr_p)
    dt_m_1bp <- expand_hits(dhs_gr, gr_m)

    res <- dt_dhs[, {
        res_p <- get_robust_edge(
            pos = dt_p_1bp[.(.BY$DHS_id), pos, nomatch = NULL],
            score = dt_p_1bp[.(.BY$DHS_id), score, nomatch = NULL],
            region_start = start, region_end = end,
            strand = "+", bg_ratio = 0.1, min_width = 20
        )
        res_m <- get_robust_edge(
            pos = dt_m_1bp[.(.BY$DHS_id), pos, nomatch = NULL],
            score = dt_m_1bp[.(.BY$DHS_id), score, nomatch = NULL],
            region_start = start, region_end = end,
            strand = "-", bg_ratio = 0.1, min_width = 20
        )
        .(
            DHS = DHS,
            center = start + floor((end - start) / 2),
            tss_plus_abs = res_p$tss,
            cov_plus = res_p$cov,
            tss_minus_abs = res_m$tss,
            cov_minus = res_m$cov
        )
    }, by = .(DHS_id)]
    res[, DHS_id := NULL]
    res[, condition := condition_name]
    res[, enhancer_width := abs(tss_plus_abs - tss_minus_abs)]
    return(res)
}

unique_conditions <- unique(meta_nascent$base_condition)

nascent_results <- mclapply(unique_conditions, function(cond) {
    cat("Analyzing pooled condition:", cond, "\n")
    sub_meta <- meta_nascent[meta_nascent$base_condition == cond, ]
    calculate_nascent_edges(
        plus_bw = sub_meta$plus_bw,
        minus_bw = sub_meta$minus_bw,
        dhs_gr = dhs_intergenic,
        condition_name = cond
    )
}, mc.cores = length(unique_conditions))

master_nascent_table <- rbindlist(nascent_results, use.names = TRUE, fill = TRUE)
fwrite(master_nascent_table, file.path(out_dir, "tss_edges_per_DHS.tsv"), sep = "\t")


#### TSS edges as bedGraph tracks (0-based)

master_nascent_table[, chrom := sub("_.*", "", DHS)]

dt_plus <- master_nascent_table[!is.na(tss_plus_abs), .(
    chrom = chrom,
    chromStart = tss_plus_abs - 1,
    chromEnd = tss_plus_abs,
    score = 1,
    condition = condition
)]

dt_minus <- master_nascent_table[!is.na(tss_minus_abs), .(
    chrom = chrom,
    chromStart = tss_minus_abs - 1,
    chromEnd = tss_minus_abs,
    score = -1,
    condition = condition
)]

dt_bedgraph <- rbindlist(list(dt_plus, dt_minus))
setorder(dt_bedgraph, condition, chrom, chromStart)

for (cond in unique(dt_bedgraph$condition)) {
    dt_sub <- dt_bedgraph[condition == cond, .(chrom, chromStart, chromEnd, score)]
    out_file <- paste0(cond, "_TSS_edges.bedGraph")
    fwrite(dt_sub, file = file.path(out_dir, out_file), sep = "\t", col.names = FALSE, quote = FALSE, scipen = 999)
    cat("Wrote", nrow(dt_sub), "TSS edges to", out_file, "\n")
}


#### Convergent -> divergent switch between the two conditions

## one row per DHS: TSSs and coverage in the control (ctrl_) and differentiated (diff_) condition
tss_shift <- master_nascent_table %>%
    dplyr::filter(condition %in% c(set$control, set$treated)) %>%
    mutate(condition = ifelse(condition == set$control, "ctrl", "diff")) %>%
    dplyr::select(DHS, tss_plus_abs, tss_minus_abs, condition, cov_plus, cov_minus) %>%
    reshape2::melt(id.vars = c("DHS", "condition")) %>%
    reshape2::dcast(DHS ~ condition + variable, value.var = "value") %>%
    drop_na() %>%
    mutate(d.ctrl = ctrl_tss_plus_abs - ctrl_tss_minus_abs,
           d.diff = diff_tss_plus_abs - diff_tss_minus_abs) %>%
    mutate(switch = d.ctrl < 0 & d.diff > 0) %>%
    dplyr::filter(ctrl_cov_plus > set$min_cov & ctrl_cov_minus > set$min_cov &
                  diff_cov_plus > set$min_cov & diff_cov_minus > set$min_cov)

fwrite(tss_shift, file.path(out_dir, "enhancer_tss_shift.tsv"), sep = "\t")

tss_shift_summary <- data.frame(
    dataset = dataset,
    control = set$control,
    differentiated = set$treated,
    min_cov = set$min_cov,
    n_enhancers = nrow(tss_shift),
    n_switch = sum(tss_shift$switch),
    percent_switch = 100 * mean(tss_shift$switch)
)
print(tss_shift_summary)
write.table(tss_shift_summary, file.path(out_dir, "enhancer_tss_shift_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

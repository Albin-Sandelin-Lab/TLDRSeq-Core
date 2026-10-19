## Enhancer DHSs: TLDR-Seq read coverage around the DHS (RBM7), ChIP-seq / chromatin signal
## per TSS distance quantile, and RBM7 transcript length vs distance between the two TSSs.
##
## Usage (from 07-Enhancer/Enhancer-analysis/): Rscript 01_enhancer_coverage_and_chip.r
##
## Enhancer DHSs with a TSS on both strands. TSS center = middle of the + and - strand TSSs;
## TSS distance (dist.tsss) = TSS+ - TSS- (negative: the two transcripts overlap / converge).
## DHSs are sorted by TSS distance and split into 9 quantiles (ntile.dist).
##
## Input : ENHANCER_ANNOTATION_RDS, DIR_READS_PER_DHS, SCALE_FACTORS_CSV (config/paths.R)
##         bigWig tracks: ENHANCER_DATA_DIR/*.bw, HISTONE_MARKS_DIR (H3K27ac, H3K4me3),
##         REMAP_BW_DIR (ReMap TFs), normalized by the area under each track (REMAP_NORM_FACTORS)
## Output (in $ENHANCER_DIR/coverage_and_chip/):
##   DHS_reads_df.tsv                     TLDR-Seq reads of the enhancer DHSs (WT, RBM7, ZCCHC8)
##   tss_distance.tsv                     TSS+, TSS-, TSS center, TSS distance, width per DHS
##   DHS_ntile_distance.tsv               TSS distance and width quantiles (9) per DHS
##   all_plus_strand_coverage.png         RBM7 + strand read coverage around the TSS center
##   all_minus_strand_coverage.png        RBM7 - strand read coverage
##   five_both_strand_coverage.png        RBM7 read 5' end (21 nt) coverage, both strands
##   metadata_summary_dt.tsv              mean signal of each track per position (+/- 1 kb from
##                                        the TSS center) and TSS distance quantile
##   metadata_summary_dt_tss.tsv          same, anchored on the + strand TSS and on the - strand
##                                        TSS (column anchor)
##   tss_center_summary.tsv               mean / median / sd TSS distance per quantile
##   STARRSeq_metaplots_Quantiles_1_2_3.pdf  median STARR-seq signal (smoothed) around the TSS
##                                        center of convergent DHSs (TSS distance < 0)
##   RBM7_split_panels_length_and_TSS_2.pdf  RBM7 mean read length per DHS and TSS distance,
##                                        per TSS distance quantile

library(rtracklayer)
library(GenomicRanges)
library(parallel)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(ggridges)
library(patchwork)

source("../../config/paths.R")
source("../../lib/tldr_reads.R")

out_dir <- file.path(ENHANCER_DIR, "coverage_and_chip")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#### bigWig tracks and normalization factors

## norm_factor = 1e8 / area under the track
norm_factor.df <- read.table(REMAP_NORM_FACTORS, sep = "\t", header = TRUE) %>%
    mutate(norm_factor = (1 / area) * 100000000) %>%
    mutate(TF = file %>% str_replace("\\.bw", "") %>% str_replace("^Pooled_", ""))

bw.files <- list.files(c(ENHANCER_DATA_DIR, REMAP_BW_DIR), pattern = "\\.bw$", full.names = TRUE) %>%
    c(file.path(HISTONE_MARKS_DIR, c("H3K27ac-mm_norm_mm39.bw", "H3K4me3-mm_norm_mm39.bw")))

names.bw <- basename(bw.files) %>% str_replace("\\.bw$", "") %>% str_replace("Pooled_", "")


#### Enhancer DHSs and their reads

enhancer.df <- readRDS(ENHANCER_ANNOTATION_RDS)

dhs.list <- enhancer.df %>%
    .$DHS

replicate.list <- c("mES_WT_rep1", "mES_WT_rep2", "mES_WT_rep3", "mES_WT_rep4", "mES_WT_rep5", "mES_WT_rep6",
                    "mES_RBM7_rep1", "mES_RBM7_rep2", "mES_RBM7_rep3",
                    "mES_ZCCHC8_rep1", "mES_ZCCHC8_rep2", "mES_ZCCHC8_rep3")

DHS_reads.list <- mclapply(dhs.list,
                           read_TLDR_reads_by_list_of_DHS,
                           query_list_replicates = replicate.list,
                           mc.cores = 50)

DHS.reads.df <- rbindlist(DHS_reads.list)
fwrite(DHS.reads.df, file.path(out_dir, "DHS_reads_df.tsv"), sep = "\t")
rm(DHS_reads.list, DHS.reads.df)
gc()

all_reads_dt <- fread(file.path(out_dir, "DHS_reads_df.tsv"))

setnames(all_reads_dt,
         c("V2", "V3", "V5", "V6", "V10", "V11", "V13", "V14"),
         c("start", "end", "pA_len", "strand", "n.exon", "exon.lengths", "DHS", "rep"))


#### TSS center and TSS distance of each DHS

tss_center.df <- enhancer.df[, c("DHS", "mean_TSS_+", "mean_TSS_-", "DHS_width")] %>%
    setNames(c("DHS", "TSS_plus", "TSS_minus", "DHS_width")) %>%
    drop_na() %>%
    mutate(tss_center = round((TSS_plus + TSS_minus) / 2)) %>%
    mutate(offset = floor(tss_center - TSS_plus)) %>%
    mutate(dist.tsss = floor(TSS_plus - TSS_minus)) %>%
    dplyr::select(DHS, tss_center, offset, dist.tsss, TSS_plus, TSS_minus, DHS_width)

write.table(tss_center.df, file.path(out_dir, "tss_distance.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## TSS distance and width quantiles
dhs.order.for.center <- tss_center.df %>%
    mutate(ntile.width = ntile(DHS_width, 9)) %>%
    mutate(ntile.dist = ntile(dist.tsss, 9)) %>%
    dplyr::select(DHS, ntile.width, ntile.dist, dist.tsss) %>%
    as.data.table

write.table(dhs.order.for.center, file.path(out_dir, "DHS_ntile_distance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


#### Reads: chrM scale factors, one row per exon

df_scale_factor <- read.csv(SCALE_FACTORS_CSV) %>%
    mutate(Rep = sample %>% str_replace("scale_", "")) %>%
    dplyr::select(scale_factor, rep = Rep) %>%
    mutate(scale = 1 / scale_factor) %>%
    dplyr::select(-scale_factor) %>%
    as.data.table()

all_reads_dt[df_scale_factor, scale := i.scale, on = "rep"]
all_reads_dt <- all_reads_dt[DHS %in% tss_center.df$DHS, ]

# Extract 'geno' (e.g., "WT", "RBM7", "ZCCHC8") from the 'rep' column (mES_WT_rep1)
all_reads_dt[, geno := tstrsplit(rep, "_")[[2]]]

# Join the reads with the tss_center dataframe (adds tss_center, offset, dist.tsss)
setkey(all_reads_dt, DHS)
tss_center_dt <- as.data.table(tss_center.df)
setkey(tss_center_dt, DHS)
all_reads_dt <- all_reads_dt[tss_center_dt, nomatch = NULL]

setnames(all_reads_dt,
         old = c("V1", "n.exon", "exon.lengths", "V12"),
         new = c("seqnames", "blockCount", "blockSizes", "blockStarts"),
         skip_absent = TRUE)

setkey(all_reads_dt, DHS, geno, strand)

# Separate unspliced reads and assign exact coordinates
unspliced <- all_reads_dt[blockCount == 1]
unspliced[, `:=`(exon_start = start, exon_end = end)]

# Isolate spliced reads
spliced <- all_reads_dt[blockCount > 1]

if (nrow(spliced) > 0) {
  # 1. Split the comma-separated strings (vectorized over the whole column)
  size_list <- strsplit(spliced$blockSizes, ",")
  start_list <- strsplit(spliced$blockStarts, ",")
  # 2. Count how many exons each read has
  n_exons <- lengths(size_list)
  # 3. Expand the data.table by duplicating rows based on exon count
  spliced_exp <- spliced[rep(1:.N, n_exons)]
  # 4. Calculate actual coordinates using fast unlisted vectors
  spliced_exp[, exon_start := start + as.numeric(unlist(start_list))]
  spliced_exp[, exon_end := exon_start + as.numeric(unlist(size_list)) - 1]
  # 5. Combine everything back into one flat table
  flat_reads_dt <- rbindlist(list(unspliced, spliced_exp), use.names = TRUE, fill = TRUE)
} else {
  flat_reads_dt <- unspliced
}


#### Read coverage per DHS, genotype and strand, relative to the TSS center

## scaled coverage of the reads (or exons) of one DHS / genotype / strand
compute_coverage <- function(reads_dt) {
    setkey(reads_dt, DHS, geno, strand)
    unique_groups <- reads_dt[, .N, keyby = .(DHS, geno, strand)][, .(DHS, geno, strand)]
    indices <- 1:nrow(unique_groups)
    results_list <- mclapply(indices, function(i) {
        ## Get the "work order" for this specific worker
        current_group <- unique_groups[i]
        ## Pull the chunk directly from shared memory (reads_dt is keyed)
        chunk <- reads_dt[current_group, nomatch = NULL]
        if (nrow(chunk) == 0) return(NULL)
        this_center <- chunk$tss_center[1]
        this_chrom  <- chunk$seqnames[1]
        gr <- GRanges(
            seqnames = this_chrom,
            ranges = IRanges(start = chunk$exon_start, end = chunk$exon_end),
            scale = chunk$scale
        )
        rle_list <- coverage(gr, weight = "scale")
        if (this_chrom %in% names(rle_list)) {
            rle_obj <- rle_list[[this_chrom]]
            res <- data.table(
                score = runValue(rle_obj),
                p_start = start(rle_obj),
                p_end = end(rle_obj)
            )[score > 0]
            if (nrow(res) > 0) {
                expanded <- res[, .(
                    position = seq(p_start, p_end),
                    score = score
                ), by = p_start]
                expanded[, relative_pos := position - this_center]
                return(data.table(
                    DHS = current_group$DHS,
                    geno = current_group$geno,
                    strand = current_group$strand,
                    position = expanded$position,
                    score = expanded$score,
                    relative_pos = expanded$relative_pos
                ))
            }
        }
        return(NULL)
    }, mc.cores = 50)
    rbindlist(results_list)
}

raw_plot_data <- compute_coverage(flat_reads_dt)
rm(flat_reads_dt)
gc()

# Define the full universe of variables to ensure no gaps in the heatmap
all_dhs_ids <- unique(tss_center.df$DHS)
all_genos <- c("WT", "RBM7", "ZCCHC8")
all_strands <- c("+", "-")

plot_data <- raw_plot_data %>%
  # Average across replicates (6 for WT, 3 for mutants)
  mutate(score = ifelse(geno == "WT", score / 6, score / 3)) %>%
  # log10 transformation; - strand coverage negative (+ up, - down)
  mutate(
    base_log_score = log10(1 + score),
    plot_score = ifelse(strand == "-", -base_log_score, base_log_score)
  ) %>%
  # Fill in missing coordinates with zeros to ensure smooth, gapless heatmaps
  tidyr::complete(DHS = all_dhs_ids,
                  geno = all_genos,
                  strand = all_strands,
                  fill = list(relative_pos = 0, plot_score = 0, score = 0, base_log_score = 0)) %>%
  drop_na()

# Only RBM7 is plotted
plot_data$geno <- factor(plot_data$geno, levels = c("RBM7"))
plot_data <- plot_data %>% drop_na()

max_dist <- 1000

## DHSs sorted by TSS distance
dhs_order <- tss_center.df %>%
  arrange(desc(dist.tsss)) %>%
  .$DHS

plot_data_filtered <- plot_data %>%
  filter(abs(relative_pos) <= max_dist) %>%
  mutate(DHS = factor(DHS, levels = dhs_order)) %>%
  filter(!is.na(DHS))

plot_data_filtered$plot_score <- ifelse(
  abs(plot_data_filtered$plot_score) > 2,
  sign(plot_data_filtered$plot_score) * 2,
  plot_data_filtered$plot_score
)

## TSS positions relative to the TSS center
tss_annotations <- tss_center.df %>%
  filter(DHS %in% dhs_order) %>%
  mutate(
    # MUST be factored exactly like the main data to align on the y-axis
    DHS = factor(DHS, levels = dhs_order),
    rel_TSS_plus = TSS_plus - tss_center,
    rel_TSS_minus = TSS_minus - tss_center
  ) %>%
  drop_na(DHS)

DHS.samples <- unique(plot_data_filtered$DHS)

plot_test <- plot_data_filtered %>% dplyr::filter(DHS %in% DHS.samples)


#### Heatmaps of the read coverage, one strand

plot_strand_heatmap <- function(plot_test, strand_choice, DHS.samples, tss_annotations, out_dir, file_prefix = "DHS") {
    rdbu <- rev(brewer.pal(n = 11, name = "RdBu"))
    ## Set dynamic titles and filenames based on the strand choice
    if (strand_choice == "+") {
        plot_title <- "+ Strand Read Coverage at DHS Sites"
        file_name <- paste0(file_prefix, "_plus_strand_coverage.png")
    } else if (strand_choice == "-") {
        plot_title <- "- Strand Read Coverage at DHS Sites"
        file_name <- paste0(file_prefix, "_minus_strand_coverage.png")
    } else {
        stop("strand_choice must be '+' or '-'")
    }
    p <- ggplot(plot_test %>% dplyr::filter(strand == strand_choice, DHS %in% DHS.samples), aes(x = relative_pos, y = DHS)) +
        geom_tile(aes(fill = plot_score)) +
        scale_fill_gradientn(
            colours = c(rdbu[2], rdbu[3], rdbu[4], "white", "white", rdbu[8], rdbu[9], rdbu[10]),
            values = scales::rescale(c(-2, -1.7, -0.6, -0.1, 0.1, 0.6, 1.7, 2)),
            breaks = c(-2, -1, 0, 1, 2),
            limits = c(-2, 2)
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
        labs(
            x = "Position relative to TSS Center (bp)",
            y = "Enhancer DHS Sites (Sorted by TSS Distance)",
            title = plot_title,
            fill = "Log10\nCoverage"
        ) +
        geom_point(data = tss_annotations %>% dplyr::filter(DHS %in% DHS.samples),
                   aes(x = rel_TSS_plus, y = DHS),
                   shape = 124, size = 0.1, color = "black") +
        geom_point(data = tss_annotations %>% dplyr::filter(DHS %in% DHS.samples),
                   aes(x = rel_TSS_minus, y = DHS),
                   shape = 124, size = 0.2, color = "#a3a3a3") +
        scale_x_continuous(labels = scales::comma, breaks = c(-1000, 0, 1000)) +
        theme_minimal() +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            panel.spacing = unit(1, "lines"),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            strip.background = element_rect(fill = "grey90", color = "grey90"),
            strip.text = element_text(face = "bold", size = 10)
        )
    save_path <- file.path(out_dir, file_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = save_path,
           plot = p,
           width = 5,
           height = 10)
    return(save_path)
}

plot_strand_heatmap(
  plot_test = plot_test,
  strand_choice = "-",
  DHS.samples = DHS.samples,
  tss_annotations = tss_annotations,
  out_dir = file.path(out_dir),
  file_prefix = "all"
)

plot_strand_heatmap(
  plot_test = plot_test,
  strand_choice = "+",
  DHS.samples = DHS.samples,
  tss_annotations = tss_annotations,
  out_dir = file.path(out_dir),
  file_prefix = "all"
)


#### Heatmap of the read 5' end coverage, both strands

## 5' end of each read: first 21 nt
flat_reads.five_dt <- copy(all_reads_dt)
flat_reads.five_dt[, blockCount := 1]
flat_reads.five_dt[, `:=`(exon_start = start, exon_end = end)]
flat_reads.five_dt[strand == "+", exon_end := exon_start + 20]
flat_reads.five_dt[strand == "-", exon_start := exon_end - 20]

raw_plot_data.five <- compute_coverage(flat_reads.five_dt)
rm(flat_reads.five_dt)
gc()

plot_data.five <- raw_plot_data.five %>%
  # Average across replicates (6 for WT, 3 for mutants)
  mutate(score = ifelse(geno == "WT", score / 6, score / 3)) %>%
  mutate(
    base_log_score = log10(1 + score),
    plot_score = ifelse(strand == "-", -base_log_score, base_log_score)
  ) %>%
  tidyr::complete(DHS = all_dhs_ids,
                  geno = all_genos,
                  strand = all_strands,
                  fill = list(position = 0, relative_pos = 0, plot_score = 0, score = 0, base_log_score = 0))

plot_data.five$geno <- factor(plot_data.five$geno, levels = c("RBM7"))
plot_data.five <- plot_data.five %>% drop_na()

plot_data_filtered.five <- plot_data.five %>%
  filter(abs(relative_pos) <= max_dist) %>%
  mutate(DHS = factor(DHS, levels = dhs_order)) %>%
  filter(!is.na(DHS))

plot_data_filtered.five$plot_score <- ifelse(
    abs(plot_data_filtered.five$plot_score) > 2,
    sign(plot_data_filtered.five$plot_score) * 2,
    plot_data_filtered.five$plot_score
)

DHS.samples.five <- unique(plot_data_filtered.five$DHS)
plot_test.five <- plot_data_filtered.five %>% dplyr::filter(DHS %in% DHS.samples.five)

## one row per DHS and strand (+ above -)
plot_both_strands_heatmap <- function(plot_test, DHS.samples, tss_annotations, out_dir, file_prefix = "DHS") {
    plot_title <- "+/- Strand Read Coverage at DHS Sites"
    file_name <- paste0(file_prefix, "_both_strand_coverage.png")
    rdbu <- rev(brewer.pal(n = 11, name = "RdBu"))
    order.dhs.strand <- paste0(matrix(rep(levels(plot_test$DHS), 2), ncol = 2, byrow = FALSE) %>% t %>% c,
                               c(".+", ".-"))
    plot_test <- plot_test %>%
        mutate(DHS.strand = paste(DHS, strand, sep = ".")) %>%
        mutate(DHS.strand = factor(DHS.strand, levels = order.dhs.strand))
    p <- ggplot(plot_test %>% dplyr::filter(DHS %in% DHS.samples), aes(x = relative_pos, y = DHS.strand)) +
        geom_tile(aes(fill = plot_score)) +
        scale_fill_gradientn(
            colours = c(rdbu[2], rdbu[3], rdbu[4], "white", "white", rdbu[8], rdbu[9], rdbu[10]),
            values = scales::rescale(c(-2, -1.7, -0.1, -0.05, 0.05, 0.1, 1.7, 2)),
            breaks = c(-2, -1, 0, 1, 2),
            limits = c(-2, 2),
            na.value = "white"
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
        labs(
            x = "Position relative to TSS Center (bp)",
            y = "Enhancer DHS Sites (Sorted by TSS Distance)",
            title = plot_title,
            fill = "Log10\nCoverage"
        ) +
        scale_x_continuous(labels = scales::comma, breaks = c(-1000, 0, 1000)) +
        theme_minimal() +
        theme(panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid = element_blank(),
              panel.spacing = unit(1, "lines"),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              strip.background = element_rect(fill = "grey90", color = "grey90"),
              strip.text = element_text(face = "bold", size = 10)
        )
    save_path <- file.path(out_dir, file_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = save_path,
           plot = p,
           width = 5,
           height = 10)
    return(save_path)
}

plot_both_strands_heatmap(
  plot_test = plot_test.five,
  DHS.samples = DHS.samples,
  tss_annotations = tss_annotations,
  out_dir = file.path(out_dir),
  file_prefix = "five"
)


#### bigWig signal around the TSS center and the two TSSs (+/- 1 kb)

## per-base signal of one track in the windows around center_df$tss_center
extract_bw_to_heatmap <- function(bw_path, center_df, window_size = 1000) {
    ## 1. Define the exact genomic windows we need
    chroms <- tstrsplit(center_df$DHS, "_")[[1]]
    gr_windows <- GRanges(
        seqnames = chroms,
        ranges = IRanges(start = center_df$tss_center - window_size,
                         end = center_df$tss_center + window_size),
        DHS = center_df$DHS,
        tss_center = center_df$tss_center
    )
    ## 2. Targeted Import: ONLY load data overlapping our windows
    bw_gr <- import.bw(bw_path, selection = BigWigSelection(gr_windows))
    ## 3. Map the extracted scores back to the specific DHS windows
    hits <- findOverlaps(bw_gr, gr_windows)
    raw_dt <- data.table(
        DHS = gr_windows$DHS[subjectHits(hits)],
        center = gr_windows$tss_center[subjectHits(hits)],
        p_start = start(bw_gr)[queryHits(hits)],
        p_end = end(bw_gr)[queryHits(hits)],
        score = bw_gr$score[queryHits(hits)]
    )
    if (nrow(raw_dt) == 0) return(NULL)
    ## 4. Trim segments that spill outside our strict window
    raw_dt[, true_start := pmax(p_start, center - window_size)]
    raw_dt[, true_end := pmin(p_end, center + window_size)]
    raw_dt <- raw_dt[true_start <= true_end]
    ## 5. Vectorized Expansion
    lens <- (raw_dt$true_end - raw_dt$true_start) + 1
    expanded_dt <- raw_dt[rep(1:.N, lens)]
    expanded_dt[, position := true_start + sequence(lens) - 1]
    expanded_dt[, relative_pos := position - center]
    ## 6. Aggregate duplicates (in case of rare overlapping BigWig bins)
    res_dt <- expanded_dt[, .(score = mean(score)), by = .(DHS, relative_pos)]
    cat(basename(bw_path))
    return(res_dt)
}

#### STARR-seq signal on convergent DHSs (TSS distance < 0)

STARRSeq_convergent.dt <- extract_bw_to_heatmap(file.path(ENHANCER_DATA_DIR, "STARRseq_merged.bw"),
                                                center_df = tss_center.df %>%
                                                    dplyr::filter(dist.tsss < 0),
                                                window_size = 1000)

## median signal per position, smoothed over 51 nt
STARRSeq_summary_dt <- STARRSeq_convergent.dt[, .(plot_score = median(score), sd.score = sd(score)), by = .(relative_pos)]

STARRSeq_summary_dt[, `:=`(plot_score_smooth = frollmean(plot_score, n = 51, align = "center", fill = NA),
                           plot_sd_smooth = frollmean(sd.score, n = 51, align = "center", fill = NA))]

g <- ggplot(STARRSeq_summary_dt %>% drop_na(), aes(x = relative_pos, y = plot_score_smooth)) +
    geom_line(size = 0.8, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
        x = "Position relative to TSS Center (bp)",
        y = "Median STARR-seq Signal (smoothed)",
        title = "Median STARR-Seq Signal on convergent DHSs (Q1, Q2, Q3)"
    ) +
    scale_x_continuous(labels = scales::comma) +
    theme_minimal() +
    theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right"
    )

ggsave(filename = file.path(out_dir, "STARRSeq_metaplots_Quantiles_1_2_3.pdf"),
       plot = g, width = 14, height = 10)

## anchored on the TSS center
list_of_dts <- mclapply(bw.files, function(bw_path) {
    extract_bw_to_heatmap(bw_path, center_df = tss_center.df, window_size = 1000)
}, mc.cores = 5)
names(list_of_dts) <- names.bw

## anchored on the + strand TSS and on the - strand TSS
list_of_dts.plus <- mclapply(bw.files, function(bw_path) {
    extract_bw_to_heatmap(bw_path, center_df = tss_center.df %>%
                                       mutate(tss_center = round(TSS_plus)),
                          window_size = 1000)
}, mc.cores = 5)

list_of_dts.minus <- mclapply(bw.files, function(bw_path) {
    extract_bw_to_heatmap(bw_path, center_df = tss_center.df %>%
                                       mutate(tss_center = round(TSS_minus)),
                          window_size = 1000)
}, mc.cores = 5)

names(list_of_dts.plus) <- names.bw
names(list_of_dts.minus) <- names.bw


#### Mean signal per position and TSS distance quantile

meta_dt <- as.data.table(dhs.order.for.center)[, .(DHS, ntile.dist)]

## number of DHSs per quantile: DHSs without signal count as 0 in the mean
quantile_counts <- meta_dt[, .(total_dhs = .N), by = ntile.dist]

norm_factor.dt <- as.data.table(norm_factor.df)[, .(mark = TF, norm_factor)]

## summed signal per quantile and position, for each track
summarize_anchors <- function(dt_list, anchor_name = NA) {
  res_list <- mclapply(names.bw, function(mark) {
    dt <- dt_list[[mark]]
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    # A. Merge the sparse data with the quantile info
    merged_dt <- dt[meta_dt, on = "DHS", nomatch = 0]
    # B. Aggregate the signal per quantile and position
    summary_dt <- merged_dt[, .(total_score = sum(score, na.rm = TRUE)),
                            by = .(ntile.dist, relative_pos = relative_pos)]
    return(summary_dt)
  }, mc.cores = 10)
  names(res_list) <- names.bw
  res_list <- res_list[!sapply(res_list, is.null)]
  res_dt <- rbindlist(res_list, idcol = "mark")
  res_dt[, anchor := anchor_name]
  return(res_dt)
}

## mean signal = summed signal / number of DHSs in the quantile, normalized per track
## (x norm_factor: every track scaled to the same total signal)
mean_signal <- function(final_summary_dt) {
    final_summary_dt <- final_summary_dt[quantile_counts, on = "ntile.dist"]
    final_summary_dt[, true_mean_score := total_score / total_dhs]
    final_summary_dt <- final_summary_dt[norm_factor.dt, on = "mark"]
    final_summary_dt[, true_mean_score := true_mean_score * norm_factor]
    final_summary_dt[, plot_score := true_mean_score]
    final_summary_dt[, ntile.dist := as.factor(ntile.dist)]
    final_summary_dt
}

tss_center.df.summarized <- tss_center.df %>%
    left_join(meta_dt, by = "DHS") %>%
    group_by(ntile.dist) %>%
    summarize(mean_tss_distance = mean(dist.tsss),
              median_tss_distance = median(dist.tsss),
              sd_tss_distance = sd(dist.tsss),
              .groups = 'drop') %>%
    mutate(ntile.dist = factor(ntile.dist))

write.table(tss_center.df.summarized, file.path(out_dir, "tss_center_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## anchored on the TSS center
final_summary_dt <- mean_signal(summarize_anchors(list_of_dts))
final_summary_dt.albin <- final_summary_dt[, .(mark, ntile.dist, relative_pos, plot_score)]
fwrite(final_summary_dt.albin, file.path(out_dir, "metadata_summary_dt.tsv"), sep = "\t")

## anchored on the + strand TSS and on the - strand TSS
final_summary_dt <- mean_signal(rbindlist(list(summarize_anchors(list_of_dts.plus, "plus"),
                                               summarize_anchors(list_of_dts.minus, "minus"))))
final_summary_dt.albin <- final_summary_dt[, .(mark, anchor, ntile.dist, relative_pos, plot_score)]
fwrite(final_summary_dt.albin, file.path(out_dir, "metadata_summary_dt_tss.tsv"), sep = "\t")


#### RBM7: mean read length per DHS vs TSS distance, per TSS distance quantile

## Filter for RBM7 and calculate genomic read length
plot_dt <- all_reads_dt[geno == "RBM7"]
plot_dt[, read_length := end - start]

## Create a lookup table for the quantiles
tss_lookup <- tss_annotations %>%
    mutate(ntile.dist = ntile(dist.tsss, 9)) %>%
    dplyr::select(DHS, ntile.dist) %>%
    as.data.table()

plot_dt <- merge(plot_dt, tss_lookup, by = "DHS", all.x = TRUE)
plot_dt <- plot_dt[!is.na(ntile.dist)]

## Weighted (chrM scale factor) mean read length for every DHS
dhs_means <- plot_dt[, .(
    mean_rl = sum(read_length * scale) / sum(scale)
), by = .(DHS, ntile.dist)]

## absolute TSS distance of every DHS
tss_dist.means <- tss_annotations %>%
    mutate(ntile.dist = ntile(dist.tsss, 9)) %>%
    group_by(ntile.dist) %>%
    mutate(mean_rl = abs(dist.tsss)) %>%
    as.data.table() %>%
    dplyr::select(DHS, ntile.dist, mean_rl) %>%
    mutate(ntile_label = factor(paste0("Q", ntile.dist), levels = paste0("Q", 9:1)))

## Q1 at the top of the plot, Q9 at the bottom
dhs_means[, ntile_label := factor(paste0("Q", ntile.dist), levels = paste0("Q", 9:1))]

dhs_means_filtered <- dhs_means[mean_rl < 3000]

p_ridges <- ggplot() +
    geom_density_ridges(
        data = dhs_means_filtered,
        aes(x = mean_rl, y = ntile_label, fill = ntile_label),
        alpha = 0.8,
        scale = 1.5,
        rel_min_height = 0.01,
        color = "grey30",
    ) +
    coord_cartesian(xlim = c(0, 1000)) +
    scale_fill_viridis_d(option = "mako", direction = -1, guide = "none") +
    theme_bw() +
    theme(
        axis.text.y = element_text(face = "bold", size = 12),
        panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor.x = element_blank()
    ) +
    labs(
        subtitle = "Mean Transcript Length per DHS",
        x = "Transcript Length (bp)",
        y = "Enhancer-TSS Distance Quantile"
    )

p_ridges.tssDist <- ggplot() +
    geom_density_ridges(
        data = tss_dist.means,
        aes(x = mean_rl, y = ntile_label, fill = ntile_label),
        alpha = 0.8,
        scale = 1.5,
        rel_min_height = 0.01,
        color = "grey30",
        ) +
    coord_cartesian(xlim = c(0, 1000)) +
    scale_fill_viridis_d(option = "mako", direction = -1, guide = "none") +
    theme_bw() +
    theme(
        axis.text.y = element_text(face = "bold", size = 12),
        panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor.x = element_blank()
    ) +
    labs(
        subtitle = "Tss distance per DHS",
        x = "absolute value of tss distance (bp)",
        y = "Enhancer-TSS Distance Quantile"
    )

final_plot <- p_ridges + p_ridges.tssDist +
    plot_layout(nrow = 2) +
    plot_annotation(
        title = "RBM7: Transcript Length vs. Enhancer-TSS Distance",
        theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
    )

ggsave(filename = file.path(out_dir, "RBM7_split_panels_length_and_TSS_2.pdf"),
       plot = final_plot, width = 12, height = 8)

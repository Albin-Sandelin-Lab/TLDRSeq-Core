## Effect of AT content on RBM7 / ZCCHC8 sensitivity, PROMPT vs first exonic vs first intronic
## reads: one polynomial model over the three read types, corrected for length.
##
## Usage (from 06-Genic/read-classes/): Rscript 06_AT_model_polynomial.r
##
## Gene reads (gene strand, unspliced) of the PROMPT DHSs, with the RefSeq select first exon /
## first intron of the gene whose TSS overlaps the DHS:
##   exon1 : 3' end before the start of the first intron
##   pmt   : 3' end in the first intron
## Reads are counted in 20 nt windows by the position of their 3' end from the gene TSS
## (first 1000 nt); DESeq2 log2FC (RBM7, ZCCHC8 vs WT) per DHS and window. AT content of the
## sequence from the TSS to the end of each window (cumulative AT count); windows on the first
## exon / intron junction are masked. First exonic reads: only DHSs whose first exonic pA-
## reads are nuclear (rank decile > 6, 04_cyto_nuc_ranks.r).
## Model, all read types together (windows with median read count > 3, log2FC < 0.1 set to 0):
##   log2FC ~ poly(log10(cumulative AT count), 2) + poly(log10(length), 2)
## The plotted log2FC is corrected for the length term (length set to its median).
##
## Input : $PROMPT_DIR/prompt_reads/DHS_all_reads_raw.tsv              (05-PROMPT/01)
##         $PROMPT_DIR/prompt_windows/at_log2FC_per_window_per_DHS.tsv (05-PROMPT/06)
##         $GENIC_DIR/TLDRSeq_nuclear_sup/cyto_nuc_ranks/tab_wide_rank_first_exon.tsv (04)
##         DHS_GENE_STRAND, REFSEQ_SELECT_GTF, SCALE_FACTORS_CSV, CHRM_DEPTH (config/paths.R)
## Output (in $GENIC_DIR/TLDRSeq_sup/AT_model_polynomial/):
##   at_content_per_window_per_DHS_gene.tsv  AT content per window, gene reads
##   model_universal_fit_<geno>_ysquared.pdf length-corrected log2FC vs cumulative AT count,
##                                           with marginal densities per read type

library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(tidyverse)
library(Biostrings)
library(BSgenome)
library(ggExtra)

source("../../config/paths.R")
source("../../lib/tldr_reads.R")
source("../../lib/deseq.R")
source("../../04-Normalization/scale_factors.R")

out_dir <- file.path(GENIC_DIR, "TLDRSeq_sup", "AT_model_polynomial")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## genome
genome <- Biostrings::readDNAStringSet(GENOME_FA, format = "fasta", use.names = TRUE)

dhs.gene.strand <- read.table(DHS_GENE_STRAND, sep = "\t", header = TRUE)


#### Gene reads of the PROMPT DHSs (unspliced)

all_reads_dt <- fread(file.path(PROMPT_DIR, "prompt_reads", "DHS_all_reads_raw.tsv"))

setnames(all_reads_dt,
         c("V2", "V3", "V5", "V6", "V10", "V11", "V13", "V14"),
         c("start", "end", "pA_len", "strand", "n.exon", "exon.lengths", "DHS", "rep"))

dhs.select <- dhs.gene.strand %>%
    dplyr::filter(gene_strand != "ambiguous") %>%
    .$DHS %>%
    unique

all_reads_dt2 <- all_reads_dt[, dhs.strand := paste0(DHS, ".", strand)]
all_reads_dt2 <- all_reads_dt2[DHS %in% dhs.select, ]
DHS.strand <- paste(dhs.gene.strand$DHS, dhs.gene.strand$gene_strand, sep = ".") %>%
    .[!grepl("ambiguous", .)]

all_gene_reads_dt <- all_reads_dt2[dhs.strand %in% DHS.strand, ]

all_gene_reads_dt <- all_gene_reads_dt[n.exon == 1, ] %>%
    .[, c("pA_len", "n.exon", "exon.lengths", "V12", "V9", "V7", "V8", "V4", "dhs.strand") := NULL]


#### First exon / first intron of the gene whose TSS overlaps the DHS (RefSeq select)

txdb <- makeTxDbFromGFF(REFSEQ_SELECT_GTF)

raw_exons_df <- exons(txdb, columns = c("GENEID", "exon_rank")) %>%
    as.data.frame()

clean_exons_df <- raw_exons_df %>%
    tidyr::unnest(cols = c(GENEID, exon_rank))

gene.start <- clean_exons_df %>%
    dplyr::filter(exon_rank == 1) %>%
    plyranges::as_granges() %>%
    plyranges::anchor_5p() %>%
    mutate(width = 1)

gene.intron1.start <- clean_exons_df %>%
    dplyr::filter(exon_rank == 1) %>%
    plyranges::as_granges() %>%
    plyranges::anchor_3p() %>%
    mutate(width = 1) %>%
    as.data.frame %>%
    dplyr::select(GENEID, intron1_start = start)

gene.intron1.end <- clean_exons_df %>%
    dplyr::filter(exon_rank == 2) %>%
    plyranges::as_granges() %>%
    plyranges::anchor_5p() %>%
    mutate(width = 1) %>%
    as.data.frame %>%
    dplyr::select(GENEID, intron1_end = start)

dhs.list <- all_gene_reads_dt$DHS %>% unique

dhs.granges <- dhs.list %>%
    data.frame(DHS = .) %>%
    separate(DHS, into = c("seqnames", "start", "end"), sep = "[_-]", remove = FALSE) %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    plyranges::as_granges()

dhs_gene_overlap <- findOverlaps(dhs.granges, gene.start, ignore.strand = TRUE)

dhs_gene_df <- data.frame(
    DHS = dhs.granges[queryHits(dhs_gene_overlap)]$DHS,
    gene_start = gene.start[subjectHits(dhs_gene_overlap)] %>% start(),
    gene_strand = gene.start[subjectHits(dhs_gene_overlap)] %>% strand() %>% as.character(),
    GENEID = gene.start[subjectHits(dhs_gene_overlap)]$GENEID
) %>%
    left_join(gene.intron1.start, by = "GENEID") %>%
    left_join(gene.intron1.end, by = "GENEID") %>%
    as.data.table

setkey(dhs_gene_df, DHS)
all_gene_reads_dt2 <- merge(all_gene_reads_dt, dhs_gene_df, by = "DHS", all.x = TRUE)

## exon1 reads: 3' end before the start of intron 1
all_gene_reads_dt3 <- all_gene_reads_dt2[(strand == "+" & end < intron1_start) | (strand == "-" & start > intron1_start)]

## pmt reads: 3' end after the start of intron 1 and before its end
all_gene_reads_dt4 <- all_gene_reads_dt2[(strand == "+" & end > intron1_start & end < intron1_end) |
                                         (strand == "-" & start < intron1_start & start > intron1_end)]


#### Read counts and AT content per window

calculate_window_counts <- function(all_reads_dt, tss_results) {
    ## 1. Prepare Data
    ## We merge the reads with the calculated TSS locations
    ## Note: We assume 'all_reads_dt' has columns: DHS, start, end, strand
    dt <- merge(as.data.table(all_reads_dt %>% dplyr::select(DHS, start, end, strand, DHS, rep)), tss_results, by = c("DHS", "strand"))
    ## 2. Identify the 3' end of the reads
    ## If strand is "+", the 3' end is the 'end' coordinate.
    ## If strand is "-", the 3' end is the 'start' coordinate.
    dt[, read_3p := ifelse(strand == "+", end, start)]
    ## 3. Calculate distance from TSS
    ## We use the read's own strand to determine the direction of "downstream"
    ## (+) Strand: Distance = Read_3p - TSS
    ## (-) Strand: Distance = TSS - Read_3p
    dt[, dist_from_tss := ifelse(strand == "+",
                                 read_3p - TSS,
                                 TSS - read_3p)]
    ## 4. Filter and Bin
    ## We keep only reads that fall into the first 1000bp (Windows 1-50)
    ## dist_from_tss < 0 means the read 3' end is "upstream" of the TSS (discarded)
    dt_filtered <- dt[dist_from_tss >= 0 & dist_from_tss < 1000]
    ## 5. Assign Window ID (Euclidean Division)
    ## 0-19bp   -> Window 1
    ## 20-39bp -> Window 2
    dt_filtered[, window_id := (dist_from_tss %/% 20) + 1]
    ## 6. Aggregate Counts
    ## Count reads per DHS per Window
    window_counts <- dt_filtered[, .(read_count = .N), by = .(DHS, window_id, rep, TSS, strand)]
    return(window_counts)
}

calculate_window_at <- function(tss_results, genome) {
    ## 0. VALIDATION
    if (!"exon1_len" %in% colnames(tss_results)) {
        stop("Error: tss_results must contain 'exon1_len' column for junction masking.")
    }
    ## 1. PREPARE DATA
    dt <- as.data.table(tss_results)
    dt[, seqnames := tstrsplit(DHS, "_", keep = 1)]
    n_wins <- 60
    win_size <- 20
    # Expand to 60 windows per gene
    expanded_dt <- dt[rep(1:.N, each = n_wins)]
    expanded_dt[, window_id := 1:n_wins, by = DHS]
    ## 2. CALCULATE COORDINATES
    # 0-based relative distances (critical for comparison with exon1_len)
    expanded_dt[, dist_start := (window_id - 1) * win_size]
    expanded_dt[, dist_end   := (window_id * win_size) - 1]
    expanded_dt[strand == "+", `:=`(
        win_start        = TSS + dist_start,
        win_end          = TSS + dist_end,
        window_all_start = TSS,
        window_all_end   = TSS + dist_end
    )]
    expanded_dt[strand == "-", `:=`(
        win_start        = TSS - dist_end,
        win_end          = TSS - dist_start,
        window_all_start = TSS - dist_end,
        window_all_end   = TSS
    )]
    ## 3. FETCH SEQUENCES & CALCULATE RAW AT
    windows_gr <- GRanges(expanded_dt$seqnames, IRanges(expanded_dt$win_start, expanded_dt$win_end), expanded_dt$strand)
    cumul_gr   <- GRanges(expanded_dt$seqnames, IRanges(expanded_dt$window_all_start, expanded_dt$window_all_end), expanded_dt$strand)
    expanded_dt[, AT_content    := as.numeric(letterFrequency(getSeq(genome, windows_gr), letters = "AT", as.prob = TRUE))]
    expanded_dt[, cumulative_AT := as.numeric(letterFrequency(getSeq(genome, cumul_gr), letters = "AT", as.prob = TRUE))]
    expanded_dt[, cumulative_AT.counts := as.numeric(letterFrequency(getSeq(genome, cumul_gr), letters = "AT", as.prob = FALSE))]
    expanded_dt[, width.cumulative := width(cumul_gr)]
    ## 4. CALCULATE METRICS
    ## A. 20nt Diff (Immediate gradient)
    expanded_dt[, diff_AT_20nt := AT_content - data.table::shift(AT_content, type = "lag"), by = DHS]
    ## B. 40nt Diff (Rolling Mean 2 vs Previous Rolling Mean 2)
    expanded_dt[, AT_40nt := frollmean(AT_content, n = 2, align = "right"), by = DHS]
    expanded_dt[, diff_AT_40nt := AT_40nt - data.table::shift(AT_40nt, n = 2, type = "lag"), by = DHS]
    ## C. 60nt Diff (Rolling Mean 3 vs Previous Rolling Mean 3)
    expanded_dt[, AT_60nt := frollmean(AT_content, n = 3, align = "right"), by = DHS]
    expanded_dt[, diff_AT_60nt := AT_60nt - data.table::shift(AT_60nt, n = 3, type = "lag"), by = DHS]
    ## 5. MASKING LOGIC (Filter out Junction artifacts)
    ## Flag any window that physically touches the junction
    expanded_dt[, is_junction := (dist_start < exon1_len) & (dist_end >= exon1_len)]
    ## Mask Cumulative AT if we are sitting on the junction
    expanded_dt[is_junction == TRUE, cumulative_AT := NA]
    expanded_dt[is_junction == TRUE, cumulative_AT.counts := NA]
    expanded_dt[is_junction == TRUE, width.cumulative := NA]
    ## Mask 20nt Diff: Invalid if Current OR Prev window is junction
    ## Span = 2 windows
    expanded_dt[, junc_in_2 := frollapply(is_junction, n = 2, FUN = max, align = "right", fill = 0), by = DHS]
    expanded_dt[junc_in_2 == 1, diff_AT_20nt := NA]
    ## Mask 40nt Diff: Invalid if Current 2 OR Prev 2 windows touched junction
    ## Span = 4 windows total
    expanded_dt[, junc_in_4 := frollapply(is_junction, n = 4, FUN = max, align = "right", fill = 0), by = DHS]
    expanded_dt[junc_in_4 == 1, diff_AT_40nt := NA]
    ## Mask 60nt Diff: Invalid if Current 3 OR Prev 3 windows touched junction
    ## Span = 6 windows total
    expanded_dt[, junc_in_6 := frollapply(is_junction, n = 6, FUN = max, align = "right", fill = 0), by = DHS]
    expanded_dt[junc_in_6 == 1, diff_AT_60nt := NA]
    ## 6. RETURN
    return(expanded_dt[, .(DHS, window_id, AT_content, cumulative_AT, diff_AT_20nt, diff_AT_40nt, diff_AT_60nt, cumulative_AT.counts, width.cumulative)])
}

tss_results <- all_gene_reads_dt3 %>%
    dplyr::select(DHS, TSS = gene_start, intron1_start, strand) %>%
    distinct() %>%
    mutate(exon1_len = ifelse(strand == "+",
                              intron1_start - TSS,       # (+) Length = End - Start
                              TSS - intron1_start)) %>%  # (-) Length = Start - End (TSS is higher coord)
    as.data.table

window_at_df <- calculate_window_at(tss_results, genome)
fwrite(window_at_df, file.path(out_dir, "at_content_per_window_per_DHS_gene.tsv"), sep = "\t")

window_counts <- calculate_window_counts(all_gene_reads_dt3, tss_results) %>%
    dplyr::filter(!grepl("ZC3H4|ZC3H18|ARS2", rep))

window_counts.pmt <- calculate_window_counts(all_gene_reads_dt4, tss_results) %>%
    dplyr::filter(!grepl("ZC3H4|ZC3H18|ARS2", rep))

window_counts.complete <- complete(window_counts,
                                   DHS,
                                   window_id,
                                   rep,
                                   fill = list(read_count = 0, TSS = NA, strand = NA)
                                   )

window_counts.pmt.complete <- complete(window_counts.pmt,
                                       DHS,
                                       window_id,
                                       rep,
                                       fill = list(read_count = 0, TSS = NA, strand = NA)
                                       )


#### DESeq2 per DHS and window: RBM7, ZCCHC8 vs WT

count.matrix <- window_counts.complete %>%
    dplyr::select(DHS, window_id, read_count, rep) %>%
    mutate(DHS.w = paste0(DHS, ".w", window_id)) %>%
    reshape2::dcast(DHS.w ~ rep, value.var = "read_count", fill = 0) %>%
    column_to_rownames("DHS.w")

count.matrix.pmt <- window_counts.pmt.complete %>%
    dplyr::select(DHS, window_id, read_count, rep) %>%
    mutate(DHS.w = paste0(DHS, ".w", window_id)) %>%
    reshape2::dcast(DHS.w ~ rep, value.var = "read_count", fill = 0) %>%
    column_to_rownames("DHS.w")

df.meta <- colnames(count.matrix) %>%
    data.frame(sample = .) %>%
    mutate(geno = sample %>% str_replace("mES_", "") %>% str_replace("_rep.", "")) %>%
    column_to_rownames("sample") %>%
    mutate(geno = factor(geno, c("WT", "RBM7", "ZCCHC8"))) %>%
    drop_na

deseq.results <- run_deseq(count.matrix, df.meta, candidate_geno = c("RBM7", "ZCCHC8"))
deseq.results.pmt <- run_deseq(count.matrix.pmt, df.meta, candidate_geno = c("RBM7", "ZCCHC8"))

at.log2FC.df <- deseq.results$df_lfc_shrunk %>%
    dplyr::select(DHS, log2FoldChange, metric, padj) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    separate(DHS, into = c("DHS", "window_id"), sep = "\\.w") %>%
    mutate(window_id = as.integer(window_id)) %>%
    left_join(window_at_df, by = c("DHS", "window_id")) %>%
    left_join(tss_results, by = c("DHS"))

at.log2FC.pmt.df <- deseq.results.pmt$df_lfc_shrunk %>%
    dplyr::select(DHS, log2FoldChange, metric, padj) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    separate(DHS, into = c("DHS", "window_id"), sep = "\\.w") %>%
    mutate(window_id = as.integer(window_id)) %>%
    left_join(window_at_df, by = c("DHS", "window_id")) %>%
    left_join(tss_results, by = c("DHS"))

window_counts_summary <- window_counts %>%
    separate(rep, into = c("garb", "geno", "repn"), sep = "_") %>%
    dplyr::select(-garb, -repn) %>%
    group_by(window_id, geno, DHS) %>%
    summarize(median_count = median(read_count))

window_counts_summary.pmt <- window_counts.pmt %>%
    separate(rep, into = c("garb", "geno", "repn"), sep = "_") %>%
    dplyr::select(-garb, -repn) %>%
    group_by(window_id, geno, DHS) %>%
    summarize(median_count = median(read_count))


#### PROMPT, first exonic and first intronic windows together

at.log2FC.prompts.df <- fread(file.path(PROMPT_DIR, "prompt_windows", "at_log2FC_per_window_per_DHS.tsv")) %>%
    mutate(geno = as.character(geno))

## nuclear / cytoplasmic first exonic pA- reads (rank deciles; rank 1 = most cytoplasmic)
tab.categories.exons <- read.table(file.path(GENIC_DIR, "TLDRSeq_nuclear_sup", "cyto_nuc_ranks", "tab_wide_rank_first_exon.tsv"), sep = "\t", header = TRUE) %>%
    dplyr::filter(label_3prime_pA == "first_exon_pA-") %>%
    dplyr::select(DHS, nqtile)

## first exonic windows: nuclear DHSs only (deciles 7-10; DHSs without decile are dropped)
at.log2FC.all.df <- bind_rows(
    at.log2FC.df %>% mutate(source = "exon1") %>%
    left_join(window_counts_summary, by = c("DHS", "geno", "window_id")) %>%
    dplyr::select(-intron1_start, -exon1_len) %>%
    mutate(median_count = ifelse(is.na(median_count), 0, median_count)),
    ##
    at.log2FC.pmt.df %>% mutate(source = "pmt") %>%
    left_join(window_counts_summary.pmt, by = c("DHS", "geno", "window_id")) %>%
    dplyr::select(-intron1_start, -exon1_len) %>%
    mutate(median_count = ifelse(is.na(median_count), 0, median_count)),
    ##
    at.log2FC.prompts.df %>% mutate(source = "prompts")
) %>%
    left_join(tab.categories.exons, by = "DHS") %>%
    mutate(nqtile = ifelse(nqtile %in% c(5, 6), NA, nqtile)) %>%
    mutate(location = ifelse(nqtile > 5, "nuclear", "cytosolic")) %>%
    mutate(location = ifelse(source != "exon1", "unknown", location)) %>%
    dplyr::filter(location != "cytosolic")


#### Polynomial model and plot, per genotype

for (geno.elt in c("RBM7", "ZCCHC8")) {

    window.infos.long <- at.log2FC.all.df %>%
        mutate(log2FC = ifelse(log2FoldChange < 0.1, 0, log2FoldChange)) %>%
        dplyr::filter(geno == geno.elt) %>%
        dplyr::filter(median_count > 3) %>%
        mutate(source = factor(source, levels = c("exon1", "prompts", "pmt")))

    window.infos.clean <- window.infos.long %>%
        dplyr::filter(!is.na(log2FC),
                      !is.na(cumulative_AT.counts),
                      !is.na(source))

    poly_mod <- lm(log2FC ~ poly(log10(cumulative_AT.counts), 2) + poly(log10(width.cumulative), 2),
                   data = window.infos.clean)

    ## remove the length term: prediction minus prediction at median length
    pred_all <- predict(poly_mod)
    pred_no_width <- predict(poly_mod, newdata = transform(window.infos.clean,
                                                           width.cumulative = median(width.cumulative)))

    window.infos.clean$log2FC_poly_corrected <- window.infos.clean$log2FC - (pred_all - pred_no_width)

    plot.df <- window.infos.clean %>%
        mutate(source = case_when(
                   source == "prompts" ~ "Prompts",
                   source == "exon1" ~ "First exonic",
                   source == "pmt" ~ "First intronic"
               )) %>%
        mutate(source = factor(source, levels = c("First exonic", "Prompts", "First intronic")))

    ## points per read type, one quadratic fit over all read types
    p <- ggplot(plot.df, aes(x = log10(cumulative_AT.counts),
                             y = log2FC_poly_corrected,
                             colour = source,
                             fill = source)) +
        geom_point(alpha = 0.1, size = 0.5) +
        geom_smooth(aes(group = 1),
                    method = "lm",
                    formula = y ~ poly(x, 2),
                    color = "black",
                    linetype = "dashed") +
        theme_bw() +
        labs(y = "Log2FC (Non-linear Length Correction)", x = "log10(AT Counts)",
             title = geno.elt) +
        theme(legend.position = "bottom") +
        xlim(c(1, 3))

    ## marginal densities per read type
    p_marginal <- ggMarginal(p,
                             type = "density",
                             groupColour = TRUE,
                             groupFill = TRUE,
                             alpha = 0.4)

    ggsave(file.path(out_dir, paste0("model_universal_fit_", geno.elt, "_ysquared.pdf")),
           plot = p_marginal,
           width = 8,
           height = 6)
}

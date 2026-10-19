## Change of H3K27ac and p300 at the enhancer DHSs between mESC and anterior neural progenitors
## (AntNPC, GSE89209), per TSS distance quantile.
##
## Usage (from 07-Enhancer/Enhancer-analysis/): Rscript 02_enhancer_chip_logFC.r
##
## Signal of each DHS = mean coverage +/- 250 bp around the TSS center, normalized by the area
## under the track (signal / area x 1e8: same total signal for every track).
## log2FC = log2(AntNPC + 1) - log2(mESC + 1), per DHS and mark.
##
## Input : $CHIP_GSE89209_DIR/*.bw (config/paths.R)
##         $ENHANCER_DIR/coverage_and_chip/tss_distance.tsv (01)
## Output (in $ENHANCER_DIR/chip_logFC/):
##   bw_areas.tsv                                     area under each track
##   logFC_per_DHS.tsv                                log2FC per DHS and mark
##   fraction_positive_logFC_by_distance_quantile.pdf fraction of DHSs with log2FC > 0 per TSS
##                                                    distance quantile

library(rtracklayer)
library(GenomicRanges)
library(parallel)
library(data.table)
library(tidyverse)

source("../../config/paths.R")

out_dir <- file.path(ENHANCER_DIR, "chip_logFC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


#### Design: H3K27ac and p300, mESC (control) and AntNPC

## file names: GSM2360930_WT_mESC_ChIP_H3K27me3.bw
master_design <- tibble(file_path = list.files(CHIP_GSE89209_DIR, pattern = "\\.bw$", full.names = TRUE)) %>%
    dplyr::filter(!grepl("input", file_path, ignore.case = TRUE)) %>%
    mutate(
        filename = basename(file_path),
        TF = sapply(strsplit(filename, "_"), function(x) gsub("\\.bw", "", x[5])),
        timepoint = sapply(strsplit(filename, "_"), `[`, 3), ## mESC or AntNPC
        annotation = "AntNPC_diff",
        GSE = "GSE89209",
        is_control = timepoint == "mESC"
    ) %>%
    dplyr::select(-filename) %>%
    dplyr::filter(TF %in% c("H3K27ac", "p300")) %>%
    mutate(
        mark = paste(TF, annotation, timepoint, GSE, sep = "_"),
        group.scale = paste(TF, GSE, sep = "_")
    )


#### Area under each track (normalization)

get_bw_area <- function(f) {
    gr <- import.bw(f)
    area <- sum(as.numeric(width(gr)) * gr$score)
    rm(gr); gc()
    return(data.frame(file_path = f, area = area, stringsAsFactors = FALSE))
}

areas <- mclapply(master_design$file_path, get_bw_area, mc.cores = 4) %>% bind_rows()
write.table(areas %>% mutate(file = basename(file_path)) %>% dplyr::select(file, area),
            file.path(out_dir, "bw_areas.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

master_design <- master_design %>%
    left_join(areas, by = "file_path") %>%
    mutate(norm_factor = (1 / area) * 100000000)


#### Signal around the TSS center of each DHS

tss_center.df.sample <- read.table(file.path(ENHANCER_DIR, "coverage_and_chip", "tss_distance.tsv"), sep = "\t", header = TRUE)

extract_bw_to_heatmap <- function(bw_path, center_df, window_size = 1000) {
    chroms <- tstrsplit(center_df$DHS, "_")[[1]]
    gr_windows <- GRanges(seqnames = chroms, ranges = IRanges(start = center_df$tss_center - window_size, end = center_df$tss_center + window_size), DHS = center_df$DHS, tss_center = center_df$tss_center)
    bw_gr <- import.bw(bw_path, selection = BigWigSelection(gr_windows))
    hits <- findOverlaps(bw_gr, gr_windows)
    raw_dt <- data.table(DHS = gr_windows$DHS[subjectHits(hits)], center = gr_windows$tss_center[subjectHits(hits)], p_start = start(bw_gr)[queryHits(hits)], p_end = end(bw_gr)[queryHits(hits)], score = bw_gr$score[queryHits(hits)])
    if (nrow(raw_dt) == 0) return(NULL)
    raw_dt[, true_start := pmax(p_start, center - window_size)]
    raw_dt[, true_end := pmin(p_end, center + window_size)]
    raw_dt <- raw_dt[true_start <= true_end]
    lens <- (raw_dt$true_end - raw_dt$true_start) + 1
    expanded_dt <- raw_dt[rep(1:.N, lens)]
    expanded_dt[, position := true_start + sequence(lens) - 1]
    expanded_dt[, relative_pos := position - center]
    res_dt <- expanded_dt[, .(score = mean(score)), by = .(DHS, relative_pos)]
    cat(basename(bw_path), "\n")
    return(res_dt)
}

ws <- 250

list_of_dts.center <- mclapply(master_design$file_path, function(bw_path) {
    extract_bw_to_heatmap(bw_path, center_df = tss_center.df.sample, window_size = (ws + 1))
}, mc.cores = 4)

names(list_of_dts.center) <- master_design$mark

## mean signal +/- window_size around the center (positions without signal count as 0)
get_center_values <- function(dt.mark, mark_name, window_size = 100) {
    if (is.null(dt.mark) || nrow(dt.mark) == 0) return(data.table(DHS = character(), center_score = numeric(), mark = character()))
    center_dt <- dt.mark[relative_pos >= -window_size & relative_pos <= window_size,
                         .(center_score = sum(score) / (1 + 2 * window_size)), by = DHS]
    center_dt[, mark := mark_name]
    return(center_dt)
}

center.values.dt <- lapply(master_design$mark, function(mark) {
    get_center_values(list_of_dts.center[[mark]], mark_name = mark, window_size = ws)
}) %>% rbindlist()


#### log2FC AntNPC vs mESC

## TSS distance quantiles
meta_dt <- tss_center.df.sample %>%
    mutate(ntile.dist = ntile(dist.tsss, 9)) %>%
    dplyr::select(DHS, ntile.dist)

center.values.df <- center.values.dt %>%
    tidyr::complete(DHS = unique(tss_center.df.sample$DHS), mark = unique(mark), fill = list(center_score = 0)) %>%
    left_join(master_design %>% dplyr::select(mark, norm_factor, group.scale, is_control), by = "mark") %>%
    mutate(
        center_score.norm = center_score * norm_factor,
        center_score.log = log2(pmax(center_score.norm, 0) + 1)
    ) %>%
    left_join(meta_dt, by = "DHS")

ctrl.values.df <- center.values.df %>%
    dplyr::filter(is_control == TRUE) %>%
    dplyr::select(DHS, group.scale, center_score.log, mark) %>%
    group_by(DHS, group.scale) %>% slice(1) %>% ungroup()

logFC.values.joined.df <- center.values.df %>%
    dplyr::filter(is_control == FALSE) %>%
    dplyr::select(-center_score, -norm_factor, -center_score.norm) %>%
    left_join(ctrl.values.df, by = c("DHS", "group.scale"), suffix = c(".treated", ".control")) %>%
    mutate(logFC = center_score.log.treated - center_score.log.control) %>%
    dplyr::select(DHS, group.scale, ntile.dist, logFC, mark = mark.treated)

fwrite(logFC.values.joined.df, file.path(out_dir, "logFC_per_DHS.tsv"), sep = "\t")


#### Fraction of DHSs with a positive log2FC per TSS distance quantile

plot.df <- logFC.values.joined.df %>%
    mutate(ntile.dist = factor(ntile.dist, levels = 1:9))

df.sub.fraction <- plot.df %>%
    group_by(group.scale, ntile.dist) %>%
    dplyr::summarize(n.pos = sum(logFC > 0, na.rm = TRUE), n.total = n(), frac = n.pos / n.total) %>%
    ungroup()

g <- ggplot(df.sub.fraction, aes(x = ntile.dist, y = frac, group = 1)) +
    geom_line() +
    stat_summary(fun = median, geom = "point", shape = 2, size = 3, color = "black") +
    facet_wrap(~group.scale, scales = "free_y", ncol = 2) +
    theme_bw() +
    labs(title = "Fraction of Enhancers with Positive LogFC by Distance Quantile",
         x = "Distance Quantile (ntile.dist)",
         y = "Fraction of Enhancers with Positive LogFC")
ggsave(filename = file.path(out_dir, "fraction_positive_logFC_by_distance_quantile.pdf"), plot = g, width = 10, height = 5)

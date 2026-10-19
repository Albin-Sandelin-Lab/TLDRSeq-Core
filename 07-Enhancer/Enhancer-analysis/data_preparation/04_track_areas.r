## Area under each bigWig track (sum of width x score), used to normalize the tracks in
## 01_enhancer_coverage_and_chip.r (signal / area x 1e8).
##
## Usage (from 07-Enhancer/Enhancer-analysis/data_preparation/): Rscript 04_track_areas.r
##
## Input : $REMAP_DIR/pooled_WT_bigWig/*.bw (03), ENHANCER_DATA_DIR/*.bw,
##         HISTONE_MARKS_DIR/{H3K27ac,H3K4me3}-mm_norm_mm39.bw (config/paths.R)
## Output: $REMAP_DIR/TF_normalization_factors.txt (file, area)

library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(tidyverse)

source("../../../config/paths.R")

bw_files <- list.files(REMAP_BW_DIR, pattern = "\\.bw$", full.names = TRUE)

histone_files <- file.path(HISTONE_MARKS_DIR, c("H3K27ac-mm_norm_mm39.bw", "H3K4me3-mm_norm_mm39.bw")) %>%
    c(list.files(ENHANCER_DATA_DIR, full.names = TRUE, pattern = ".bw$"))

get_bw_area <- function(file_path) {
  gr <- import.bw(file_path)
  # as.numeric() prevents integer overflow
  total_area <- sum(as.numeric(width(gr)) * gr$score)
  rm(gr)
  gc()
  return(data.frame(
    file = basename(file_path),
    area = total_area,
    stringsAsFactors = FALSE
  ))
}

num_cores <- 10
results_list <- mclapply(bw_files, FUN = get_bw_area, mc.cores = num_cores)
results_list_histones <- mclapply(histone_files, FUN = get_bw_area, mc.cores = num_cores)

norm_df <- do.call(rbind, results_list) %>%
    rbind(do.call(rbind, results_list_histones))

head(norm_df)
write.table(norm_df, REMAP_NORM_FACTORS, sep = "\t", row.names = FALSE, quote = FALSE)

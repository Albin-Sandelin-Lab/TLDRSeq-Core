#!/usr/bin/env Rscript
# starrseq_per_dhs.R -- per-DHS mean STARR-seq signal from a merged BigWig
#
# Input track: public mESC STARR-seq (Peng et al. 2020, doi:10.1186/s13059-020-02156-3; GEO GSM4261631 + GSM4261632, series GSE143544), mm9 lifted over to mm39, and replicate-mean of CPM-normalised tracks. Not redistributed with this code.

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

# env_path() reads a path from an environment variable and stops when it is unset.
source("../../02-DHS/config_and_libraries.R")

# === Paths ===================================================================
# Every path is named by an environment variable, so no machine-specific location is fixed here.
# STARRSEQ_BW (input): the merged BigWig track, mm39 coordinates, whose score column is the CPM-normalised replicate-mean STARR-seq signal.
# DHS_PROPS_FILE (input): RDS holding a data.frame with at least a character DHS column, ids of the form "chrN_start-end".
# STARRSEQ_SIGNAL_FILE (output): RDS holding a data.frame with columns DHS (ids rewritten as "chrN:start-end") and STARRSEQ_mean_signal (numeric, mean score of the track over each DHS interval, in CPM units).
STARRSEQ_BW <- env_path("STARRSEQ_BW")
DHS_PROPS_FILE <- env_path("STARRSEQ_DHS_PROPS_FILE")
OUT_FILE <- env_path("STARRSEQ_SIGNAL_FILE")

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)

# Validate paths
stopifnot(
  "STARRseq bigWig does not exist" = file.exists(STARRSEQ_BW),
  "DHS properties file does not exist" = file.exists(DHS_PROPS_FILE)
)

# === 1. Build DHS intervals ==================================================
cat("Loading DHS ids from properties file...\n")
df_props <- readRDS(DHS_PROPS_FILE)
dhs_ids <- as.character(df_props$DHS)
cat(sprintf("  DHS ids: %d\n", length(dhs_ids)))

# DHS ids look like "chr10_100465175-100466419"
parts <- str_match(dhs_ids, "^(chr[^_]+)_(\\d+)-(\\d+)$")
stopifnot(
  "Cannot parse DHS id coordinates" = all(!is.na(parts[, 1]))
)

df_dhs <- data.frame(
  DHS = dhs_ids,
  chr = parts[, 2],
  dhs_start = as.numeric(parts[, 3]),
  dhs_end = as.numeric(parts[, 4]),
  stringsAsFactors = FALSE
)

# === 2. Compute mean signal per DHS (per chromosome) =========================
message("Computing mean STARRseq signal over each DHS interval...")
chr_levels <- unique(df_dhs$chr)

starrseq_mean <- numeric(nrow(df_dhs))
for (ch in chr_levels) {
  idx <- which(df_dhs$chr == ch)
  gr <- GRanges(ch, IRanges(df_dhs$dhs_start[idx], df_dhs$dhs_end[idx]))
  vals <- import.bw(STARRSEQ_BW, which = gr, as = "NumericList")
  starrseq_mean[idx] <- vapply(vals, mean, numeric(1), na.rm = TRUE)
  cat(sprintf("  %s: %d DHS done\n", ch, length(idx)))
}

# === 3. Assemble + save ======================================================
df_starrseq <- data.frame(
  DHS = paste0(df_dhs$chr, ":", df_dhs$dhs_start, "-", df_dhs$dhs_end),
  STARRseq_mean_signal = starrseq_mean,
  stringsAsFactors = FALSE
)

cat(sprintf("  STARRseq signal: %d DHS regions\n", nrow(df_starrseq)))
cat(sprintf(
  "  Signal range: [%.4f, %.4f], median=%.4f, N_nonNA=%d\n",
  min(df_starrseq$STARRseq_mean_signal, na.rm = TRUE),
  max(df_starrseq$STARRseq_mean_signal, na.rm = TRUE),
  median(df_starrseq$STARRseq_mean_signal, na.rm = TRUE),
  sum(!is.na(df_starrseq$STARRseq_mean_signal))
))

saveRDS(df_starrseq, OUT_FILE)
cat(sprintf("  Saved: %s\n", OUT_FILE))

message("Done.")

# tf_density.R -- TF peak density over DHS on the TSS-space denominator
# density = n_hit / TSS_space_eroded, the inter-TSS interval padded by 100 bp per side.

# ── 0. Packages ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
})

# env_path() reads a path from an environment variable and stops when it is unset.
source("../../02-DHS/config_and_libraries.R")

# ── 1. Paths ─────────────────────────────────────────────────────────────────
# Every input and output is named by an environment variable, so no location is fixed here.
# DHS_FILE (input): RDS holding a data.frame, one row per enhancer DHS, with a character DHS column ("chrN_start-end") and the numeric bp columns mean_TSS_+ and mean_TSS_- (mean TSS position of the plus- and minus-strand peaks called within that DHS, mm39).
# BW_DIR (input): directory holding one "<TF>.bw" BigWig file per TF, carrying ReMap peak blocks with score > 0.
# CATEGORIES_FILE (input): TSV with a header row and the columns tf_name, category and subcategory.
# TSS_SPACE_FILE (output): RDS data.frame, one row per DHS kept, columns DHS, chr, dhs_start, dhs_end, tss_space_start, tss_space_end and tss_space_erode (bp).
# TSS_SPACE_NHIT_FILE (output): RDS integer matrix, columns = TF, values = number of peak blocks overlapping that DHS; rows are the DHS of TSS_SPACE_FILE that have a positive tss_space_erode, in that file's order, so a DHS with a non-positive width appears in the table but not in the matrix.
# TSS_SPACE_DENSITY_FILE (output): RDS numeric matrix, same rows and dimnames as the n_hit matrix, values = peak blocks per bp of tss_space_erode.
# TSS_SPACE_STATS_FILE (output): TSV, one row per TF, columns tf_name, mean_density_tss_space, n_tss_space_hit, category and subcategory.
DHS_FILE <- env_path("TF_DENSITY_DHS_FILE")
BW_DIR <- env_path("TF_DENSITY_BW_DIR")
CATEGORIES_FILE <- env_path("TF_DENSITY_CATEGORIES_FILE")
TSS_SPACE_FILE <- env_path("TF_DENSITY_TSS_SPACE_FILE")
TSS_SPACE_NHIT_FILE <- env_path("TF_DENSITY_TSS_SPACE_NHIT_FILE")
TSS_SPACE_DENSITY_FILE <- env_path("TF_DENSITY_TSS_SPACE_DENSITY_FILE")
TSS_SPACE_STATS_FILE <- env_path("TF_DENSITY_TSS_SPACE_STATS_FILE")

source("remap_tf_matrix.R")

for (out_dir in unique(dirname(c(
  TSS_SPACE_FILE,
  TSS_SPACE_NHIT_FILE,
  TSS_SPACE_DENSITY_FILE,
  TSS_SPACE_STATS_FILE
)))) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
}

# ── 2. Load data ─────────────────────────────────────────────────────────────
message("Loading enhancer DHS...")
df_enhancer_DHS <- readRDS(DHS_FILE)
message(sprintf("  %d enhancer DHS loaded.", nrow(df_enhancer_DHS)))

# Validate required columns
required_cols <- c("DHS", "mean_TSS_+", "mean_TSS_-")
missing_cols <- setdiff(required_cols, names(df_enhancer_DHS))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "df_enhancer_DHS missing required columns: %s",
    paste(missing_cols, collapse = ", ")
  ))
}

df_categories <- load_tf_categories(CATEGORIES_FILE)

# Discover TF names from the BigWig directory
all_tfs <- tools::file_path_sans_ext(
  list.files(BW_DIR, pattern = "\\.bw$", full.names = FALSE)
)
if (length(all_tfs) == 0) {
  stop(sprintf("No .bw files found in BigWig directory: %s", BW_DIR))
}
message(sprintf("  Found %d TF BigWig files in %s", length(all_tfs), BW_DIR))

# DHS id format: "chr10_100465175-100466419"
dhs_pattern <- "(chr[0-9XYMmt]+)_(\\d+)-(\\d+)"

# ── 3. Build TSS_space_eroded GRanges ────────────────────────────────────────
# TSS_space_eroded = [min(mean_TSS_+, mean_TSS_-) - 100, max(mean_TSS_+, mean_TSS_-) + 100]
# Padding the inter-TSS interval by 100 bp on each side guarantees a minimum width of 200 bp, so the density cannot blow up when the peaks coincide.
message("\nPhase 1: density_tss_space — building TSS_space region windows")

df_tss_space <- df_enhancer_DHS %>%
  select(DHS, `mean_TSS_+`, `mean_TSS_-`) %>%
  filter(!is.na(`mean_TSS_+`) & !is.na(`mean_TSS_-`))

# Parse chromosome and coordinates from the DHS string
df_tss_space <- df_tss_space %>%
  mutate(
    chr = str_extract(DHS, dhs_pattern, group = 1),
    dhs_start = as.integer(str_extract(DHS, dhs_pattern, group = 2)),
    dhs_end = as.integer(str_extract(DHS, dhs_pattern, group = 3)),
    # TSS_space_eroded: pad the inter-TSS interval by 100 bp on each side
    tss_space_start = pmin(`mean_TSS_+`, `mean_TSS_-`) - 100L,
    tss_space_end = pmax(`mean_TSS_+`, `mean_TSS_-`) + 100L,
    tss_space_erode = tss_space_end - tss_space_start
  )

df_tss_space %>% saveRDS(TSS_SPACE_FILE)

# Filter out DHS where TSS_space_eroded is zero or negative (should not happen with padding, but a missing or malformed TSS pair must not pass silently).
n_before_ts <- nrow(df_tss_space)
df_tss_space <- df_tss_space %>% filter(tss_space_erode > 0)
n_dropped_ts <- n_before_ts - nrow(df_tss_space)
if (n_dropped_ts > 0) {
  message(sprintf(
    "  Dropped %d DHS with non-positive TSS_space_eroded (missing TSS peaks).",
    n_dropped_ts
  ))
}
message(sprintf(
  "  %d DHS with valid TSS_space_eroded regions.",
  nrow(df_tss_space)
))

gr_tss_space <- GRanges(
  seqnames = df_tss_space$chr,
  ranges = IRanges(
    start = df_tss_space$tss_space_start,
    end = df_tss_space$tss_space_end
  ),
  DHS = df_tss_space$DHS
)

# ── 4. Count TF overlaps in TSS_space_eroded regions ─────────────────────────
message(
  "\nPhase 2: Counting TF overlaps in TSS_space_eroded regions..."
)

t_start_ts <- proc.time()

df_tss_space_nhit <- count_remap_tf_overlap(
  dhs_input = gr_tss_space,
  tf_names = all_tfs,
  bw_dir = BW_DIR,
  verbose = TRUE
)

elapsed_ts <- proc.time() - t_start_ts
message(sprintf("Done. Elapsed: %.1f min", elapsed_ts["elapsed"] / 60))

# ── 5. Build TSS_space_eroded density matrix ─────────────────────────────────
message("\nPhase 3: Building TSS_space_eroded density matrices...")

mat_tss_space_nhit <- as.matrix(df_tss_space_nhit[, -1])
rownames(mat_tss_space_nhit) <- df_tss_space_nhit$DHS

tss_space_widths <- df_tss_space$tss_space_erode[match(
  rownames(mat_tss_space_nhit),
  df_tss_space$DHS
)]
mat_density_tss_space <- mat_tss_space_nhit / tss_space_widths

saveRDS(
  mat_tss_space_nhit,
  TSS_SPACE_NHIT_FILE
)
saveRDS(mat_density_tss_space, TSS_SPACE_DENSITY_FILE)
message("  TSS_space density matrices saved.")

# ── 6. Summary statistics ────────────────────────────────────────────────────
message("\nPhase 4: Summary statistics...")

stats_density_tss_space <- data.frame(
  tf_name = colnames(mat_density_tss_space),
  mean_density_tss_space = colMeans(mat_density_tss_space),
  n_tss_space_hit = colSums(mat_tss_space_nhit > 0),
  stringsAsFactors = FALSE
) %>%
  left_join(df_categories, by = "tf_name") %>%
  arrange(desc(mean_density_tss_space))

write.table(
  stats_density_tss_space,
  TSS_SPACE_STATS_FILE,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n=== Per-TF density_tss_space_eroded (top 20) ===\n")
print(
  head(
    stats_density_tss_space[, c(
      "tf_name",
      "category",
      "mean_density_tss_space",
      "n_tss_space_hit"
    )],
    20
  ),
  row.names = FALSE
)

# DHS excluded from the sub-region analysis (missing TSS peaks)
n_excluded_tss_space <- nrow(df_enhancer_DHS) - nrow(df_tss_space)
cat(sprintf(
  "\n=== Sub-region coverage ===\n  DHS with TSS_space_eroded region: %d / %d (%.1f%%)  [excluded: %d]\n",
  nrow(df_tss_space),
  nrow(df_enhancer_DHS),
  100 * nrow(df_tss_space) / nrow(df_enhancer_DHS),
  n_excluded_tss_space
))

# ── 7. Session summary ───────────────────────────────────────────────────────
cat("\n=== Run complete ===\n")
cat(sprintf(
  "  DHS (density_tss_space_eroded): %d\n",
  nrow(mat_density_tss_space)
))
cat(sprintf("  TFs analysed           : %d\n", ncol(mat_density_tss_space)))
cat("  Files written:\n")
for (f in c(
  TSS_SPACE_FILE,
  TSS_SPACE_DENSITY_FILE,
  TSS_SPACE_NHIT_FILE,
  TSS_SPACE_STATS_FILE
)) {
  cat(sprintf("    %s\n", f))
}

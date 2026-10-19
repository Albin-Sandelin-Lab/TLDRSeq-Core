# 02_build_freq_table.R -- build the genic category frequency table (GENIC_FREQ_RDS)
#
# The output is the heatmap data table: one row per DHS x feature level x genotype,
# carrying feature_count, geno_count and scaled_freq for the genic heatmap
# (06-Genic/read-classes/03_read_categories_heatmap.r).
#
# The steps, in order:
#   1. Per-replicate scale factors from chrM depth (04-Normalization/scale_factors.R).
#   2. Per-read 3' labels to feature levels, each count divided by its replicate's scale
#      factor, replicates pooled within a DHS x genotype, each feature expressed as a
#      fraction of its group's pooled total (01_genic_proportion_curation.R).
#   3. Save that pooled table.
#
# Run from inside this folder:
#   Rscript 02_build_freq_table.R

if (!file.exists("00_config.R")) {
  stop(
    "Run this script from inside 06-Genic/proportions/, so that the source() calls and ",
    "the paths in 00_config.R resolve: Rscript 02_build_freq_table.R"
  )
}

source("00_config.R")
source("../../04-Normalization/scale_factors.R")
source("01_genic_proportion_curation.R")


# ---- Step 1: scale factors ------------------------------------------------
df_scale_factor <- define_scale_factor(file_depth_chrM_replicate)


# ---- Step 2: curate, rescale, pool ---------------------------------------
df_category <- readRDS(file_category_table)

list_freq <- calculate_pooled_counts_then_proportion(
  df_category,
  df_scale_factor = df_scale_factor
)


# ---- Step 3: save the pooled table ---------------------------------------
df_freq_long <- list_freq$df_freq_long
dir.create(dirname(file_genic_freq_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(df_freq_long, file_genic_freq_rds)

message(sprintf(
  "Wrote %s: %s rows, %d DHS, %d feature levels, genotypes %s",
  file_genic_freq_rds,
  format(nrow(df_freq_long), big.mark = ","),
  n_distinct(df_freq_long$DHS),
  n_distinct(df_freq_long$feature_id),
  paste(levels(df_freq_long$geno), collapse = "/")
))

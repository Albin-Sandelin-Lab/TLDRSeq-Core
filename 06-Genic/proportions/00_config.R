# 00_config.R -- inputs, output and sample universe for the heatmap data table
#
# This folder builds one table (GENIC_FREQ_RDS, config/paths.R), the input of the genic
# heatmap (06-Genic/read-classes/03_read_categories_heatmap.r): one row per
# DHS x feature level x genotype, carrying the feature's scale-factor-normalised share of
# its group's pooled reads. Every path and constant the build needs is declared here, so
# no function in this folder hides a path or a sample list inside itself.
#
# The build is three steps, run in this order by 02_build_freq_table.R:
#   1. 04-Normalization/scale_factors.R -- one scale factor per replicate, from chrM depth.
#   2. 01_genic_proportion_curation.R -- read labels to feature levels, divide by scale
#      factor, pool replicates within a DHS x genotype, express each feature as a
#      fraction of the pooled total.
#   3. 02_build_freq_table.R       -- calls 1 and 2 and saves the table.

suppressPackageStartupMessages(library(tidyverse))


# ---- Inputs ----------------------------------------------------------------
# Paths are set in config/paths.R.
source("../../config/paths.R")

# Per-read 3' classification table: one row per DHS x sample x label combination, with
# the read count n and the gene-relative 3' labels (label_tts, label_cdsts,
# label_premature_termination, label_premature, unspliced, polyA_tail,
# label_tts_falling_region) plus mean_TTS and typical_cdsss.
# Produced by 06-Genic/read-classes/label_reads_TLDRSeq_sup.R.
file_category_table <- TTS_PROP_CATEGORIES_RDS

# Per-replicate mitochondrial-chromosome depth: one line per BED12 input file, two
# whitespace-separated columns, "<depth> <path/to/bed12_file>". The replicate name is
# parsed out of the path, so the paths must carry names of the form mES_<GENOTYPE>_rep<n>.
file_depth_chrM_replicate <- CHRM_DEPTH


# ---- Output ----------------------------------------------------------------
# Frequency of each category per DHS and genotype (read by
# 06-Genic/read-classes/03_read_categories_heatmap.r).
file_genic_freq_rds <- GENIC_FREQ_RDS


# ---- Sample universe -------------------------------------------------------
# The 21 sequenced replicates, in the order the scale factors are computed and the
# sample levels are set in. Six WT replicates, three each of the five knockdowns.
# The ZC3H4 replicates carry scale factors like every other replicate, but their
# reads are dropped when the counts are curated (see build_gene_count_long() in
# 01_genic_proportion_curation.R), so they do not appear in the output table.
all_replicates <- c(
  "mES_WT_rep1", "mES_WT_rep2", "mES_WT_rep3",
  "mES_WT_rep4", "mES_WT_rep5", "mES_WT_rep6",
  "mES_ARS2_rep1", "mES_ARS2_rep2", "mES_ARS2_rep3",
  "mES_ZC3H18_rep1", "mES_ZC3H18_rep2", "mES_ZC3H18_rep3",
  "mES_ZC3H4_rep1", "mES_ZC3H4_rep2", "mES_ZC3H4_rep3",
  "mES_RBM7_rep1", "mES_RBM7_rep2", "mES_RBM7_rep3",
  "mES_ZCCHC8_rep1", "mES_ZCCHC8_rep2", "mES_ZCCHC8_rep3"
)

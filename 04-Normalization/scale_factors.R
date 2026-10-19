# scale_factors.R -- per-replicate scale factors from mitochondrial depth,
# and the replicate-name parser the counting code needs.
#
# Read counts are made comparable across replicates by dividing each replicate's counts
# by its own scale factor, and the factor is derived here from the depth of the
# mitochondrial chromosome (chrM). Where that depth table comes from is not this file's
# concern: it is passed in as a path (CHRM_DEPTH, config/paths.R).
#
# The derivation is three steps:
#   1. Read the per-replicate chrM depth table: one depth per BED12 input file, with the
#      replicate and genotype names parsed out of that file's path.
#   2. Take the mean chrM depth across all replicates. That mean is the common
#      denominator, so the factors average 1.
#   3. Divide each replicate's own chrM depth by that mean.
#
# Because the factor is a ratio of depths, a deeply sequenced replicate receives a factor
# above 1 and is therefore scaled down when its counts are divided by it.
#
# Used by 01_write_scale_factors.R (table SCALE_FACTORS_CSV), lib/deseq.R (run_deseq()) and
# 06-Genic/proportions/.


# ---- define_scale_factor ---------------------------------------------------
# Derive one scale factor per replicate from chrM depth.
#
# Input format: the chrM depth table is read with read.table() and holds two
# whitespace-separated columns per line, "<depth> <path/to/bed12_file>". The sample id is
# parsed out of the second column, so the paths must contain names of the form
# mES_<GENOTYPE>_rep<n>.
#
# Every parsed sample id is checked against the sample universe before the factors are
# built: a silently NA sample, or a sample outside the expected set, would drop out of the
# row means and quietly change every other replicate's factor.
#
# Parameters:
#   file_depth_chrM_replicate -- character. Table of per-replicate chrM depth: two
#                                whitespace-separated columns, "<depth> <path/to/bed12_file>".
#   sample_levels             -- character vector giving the expected replicates and their
#                                order (default: all_replicates, the 21 TLDRSeq_sup replicates,
#                                to be defined by the caller).
#
# Returns:
#   data.frame with columns:
#     sample        -- replicate name (e.g. "mES_WT_rep1")
#     scale_factor  -- chrM depth divided by the mean chrM depth across replicates
define_scale_factor <- function(
  file_depth_chrM_replicate,
  sample_levels = all_replicates
) {
  if (!file.exists(file_depth_chrM_replicate)) {
    stop("chrM depth file not found: ", file_depth_chrM_replicate)
  }

  df_depth_mitochondrial <- read.table(file_depth_chrM_replicate) %>%
    transmute(
      depth = V1,
      file = V2,
      sample = sapply(file, function(x) {
        str_extract(x, "(mES_\\w+_rep\\d+)_", 1)
      }),
      type = "MT",
      geno = sapply(sample, function(x) {
        str_extract(x, "mES_(\\w+)_rep", 1)
      })
    )

  # Validate every sample id parsed. A silently NA sample would drop out of the row means
  # below and quietly change every other replicate's factor.
  if (anyNA(df_depth_mitochondrial$sample) ||
        anyNA(df_depth_mitochondrial$geno)) {
    stop(
      "Could not parse a sample or genotype id from the chrM depth file. ",
      "Expected paths containing names like mES_WT_rep1. Unparsed rows: ",
      sum(is.na(df_depth_mitochondrial$sample) |
            is.na(df_depth_mitochondrial$geno))
    )
  }

  unexpected_samples <- setdiff(df_depth_mitochondrial$sample, sample_levels)
  if (length(unexpected_samples) > 0) {
    stop(sprintf(
      "chrM depth file contains sample(s) outside the expected universe: %s",
      paste(unexpected_samples, collapse = ", ")
    ))
  }

  df_depth_mito_wide <- df_depth_mitochondrial %>%
    mutate(sample = factor(sample, levels = sample_levels)) %>%
    pivot_wider(
      id_cols = -c("type", "geno", "file"),
      names_from = sample,
      values_from = depth
    )

  tmp_avg_mito_counts <- rowMeans(df_depth_mito_wide)

  df_scale_factors <- df_depth_mito_wide %>%
    transmute(
      scale_mES_WT_rep1 = mES_WT_rep1 / tmp_avg_mito_counts,
      scale_mES_WT_rep2 = mES_WT_rep2 / tmp_avg_mito_counts,
      scale_mES_WT_rep3 = mES_WT_rep3 / tmp_avg_mito_counts,
      scale_mES_WT_rep4 = mES_WT_rep4 / tmp_avg_mito_counts,
      scale_mES_WT_rep5 = mES_WT_rep5 / tmp_avg_mito_counts,
      scale_mES_WT_rep6 = mES_WT_rep6 / tmp_avg_mito_counts,

      scale_mES_ARS2_rep1 = mES_ARS2_rep1 / tmp_avg_mito_counts,
      scale_mES_ARS2_rep2 = mES_ARS2_rep2 / tmp_avg_mito_counts,
      scale_mES_ARS2_rep3 = mES_ARS2_rep3 / tmp_avg_mito_counts,

      scale_mES_RBM7_rep1 = mES_RBM7_rep1 / tmp_avg_mito_counts,
      scale_mES_RBM7_rep2 = mES_RBM7_rep2 / tmp_avg_mito_counts,
      scale_mES_RBM7_rep3 = mES_RBM7_rep3 / tmp_avg_mito_counts,

      scale_mES_ZCCHC8_rep1 = mES_ZCCHC8_rep1 / tmp_avg_mito_counts,
      scale_mES_ZCCHC8_rep2 = mES_ZCCHC8_rep2 / tmp_avg_mito_counts,
      scale_mES_ZCCHC8_rep3 = mES_ZCCHC8_rep3 / tmp_avg_mito_counts,

      scale_mES_ZC3H18_rep1 = mES_ZC3H18_rep1 / tmp_avg_mito_counts,
      scale_mES_ZC3H18_rep2 = mES_ZC3H18_rep2 / tmp_avg_mito_counts,
      scale_mES_ZC3H18_rep3 = mES_ZC3H18_rep3 / tmp_avg_mito_counts,

      scale_mES_ZC3H4_rep1 = mES_ZC3H4_rep1 / tmp_avg_mito_counts,
      scale_mES_ZC3H4_rep2 = mES_ZC3H4_rep2 / tmp_avg_mito_counts,
      scale_mES_ZC3H4_rep3 = mES_ZC3H4_rep3 / tmp_avg_mito_counts
    )

  df_scale_factors_long <- df_scale_factors %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "scale_factor"
    ) %>%
    transmute(
      sample = str_extract(sample, "mES_([A-Z0-9]+)_(rep\\d)"),
      scale_factor
    )

  return(df_scale_factors_long)
}


# ---- validate_scale_factor_table -------------------------------------------
# Assert that a scale factor table is complete and usable before it is joined onto counts.
# Every check corresponds to a way a silently wrong scaling could otherwise pass unnoticed:
# a duplicated sample picking the wrong factor, an NA factor propagating NA through the
# sum, or a non-positive factor reversing the sample's relative depth.
#
# Parameters:
#   df_scale_factor -- data.frame with columns `sample` and `scale_factor`
#
# Returns:
#   NULL invisibly; stops on the first violation.
validate_scale_factor_table <- function(df_scale_factor) {
  stopifnot(
    "df_scale_factor must be a data.frame" = is.data.frame(df_scale_factor)
  )

  missing_cols <- setdiff(c("sample", "scale_factor"), colnames(df_scale_factor))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "df_scale_factor is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (anyDuplicated(df_scale_factor$sample) > 0) {
    duplicated_samples <- unique(df_scale_factor$sample[duplicated(
      df_scale_factor$sample
    )])
    stop(sprintf(
      "define_scale_factor() returned duplicated samples: %s. Each sample must have exactly one scale factor.",
      paste(duplicated_samples, collapse = ", ")
    ))
  }

  if (anyNA(df_scale_factor$sample) || anyNA(df_scale_factor$scale_factor)) {
    stop("define_scale_factor() returned NA sample names or NA scale factors.")
  }

  if (
    any(!is.finite(df_scale_factor$scale_factor)) ||
      any(df_scale_factor$scale_factor <= 0)
  ) {
    stop(
      "define_scale_factor() returned non-finite or non-positive scale factors. Scaling requires strictly positive finite values."
    )
  }
}


# ---- extract_genotype_from_sample ------------------------------------------
# Parse the genotype out of a replicate name.
#
# Parameters:
#   sample -- character or factor, names like "mES_WT_rep1"
#
# Returns:
#   character vector of genotypes ("WT", "ARS2", ...); stops if any name does not match
#   the expected pattern.
extract_genotype_from_sample <- function(sample) {
  stopifnot(
    "sample must be character or factor" = is.character(sample) ||
      is.factor(sample)
  )

  sample <- as.character(sample)
  geno <- stringr::str_match(sample, "^mES_([A-Z0-9]+)_rep\\d+$")[, 2]

  if (anyNA(geno)) {
    bad_samples <- unique(sample[is.na(geno)])
    stop(sprintf(
      "Failed to extract genotype from sample names. Expected names like mES_WT_rep1. Problematic samples: %s",
      paste(bad_samples, collapse = ", ")
    ))
  }

  return(geno)
}

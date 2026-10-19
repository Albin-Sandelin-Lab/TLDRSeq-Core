# 01_genic_proportion_curation.R -- curate genic read classes, then express them as
# proportions of pooled scale-factor-normalised counts.
#
# This file concerns genic reads only. A genic read is one assigned to a gene and
# classified by where its 3' end falls relative to that gene's structure; the feature
# levels below are those gene-relative classes. Enhancer reads are not classified this way.
#
# Curation and proportioning are kept apart because they answer different questions.
#   1. Curation turns each read's 3' labels into one feature level plus a pA+/- suffix
#      recording whether the read carried a polyA tail.
#   2. Proportioning collapses the curated reads into one row per DHS x sample x feature
#      level, divides each count by its replicate's scale factor, pools replicates within a
#      DHS x genotype, and expresses each feature as a fraction of its group's pooled total.
#
# Requires 04-Normalization/scale_factors.R for extract_genotype_from_sample() and
# validate_scale_factor_table().


# ---- Feature levels --------------------------------------------------------
# The read-level 3' classes are produced upstream as labels (label_tts, label_cdsts,
# label_premature_termination, label_premature, unspliced, polyA_tail,
# label_tts_falling_region). This file turns those labels into the feature levels the
# proportions are taken over, and the level list below is the contract: it fixes both the
# membership and the order of the denominator.
#
# "others" is the catch-all for a read whose label combination matches no named level.
# Every read contributes to exactly one level.
SELECTED_FEATURE_LEVELS <- c(
  "single_exon_pA+",
  "single_exon_pA-",
  "first_exon_pA+",
  "first_exon_pA-",
  "gradual_degradation_pA-",
  "gradual_degradation_pA+",
  "premature_termination_pA-",
  "premature_termination_pA+",
  "UTR3_pA-",
  "UTR3_pA+",
  "over_elongation_pA-",
  "over_elongation_pA+",
  "typical_TTS_pA-",
  "others",
  "typical_TTS_pA+"
)


# ---- validate_proportion_input ---------------------------------------------
# Assert that the incoming category table carries every label column the curation reads,
# and that the counts and flags are usable: non-negative non-NA numeric counts, logical
# non-NA flags, no NA in any label column, and strand only "+" or "-".
#
# An NA label would otherwise fall through every case_when() branch and land silently in
# "others", changing the denominator without any error.
#
# Parameters:
#   df_category -- data.frame, the per-read 3' classification table
#
# Returns:
#   NULL invisibly; stops on the first violation.
validate_proportion_input <- function(df_category) {
  required_cols <- c(
    "DHS",
    "sample",
    "strand",
    "mean_TTS",
    "typical_cdsss",
    "label_tts_falling_region",
    "label_premature_termination",
    "label_tts",
    "label_premature",
    "unspliced",
    "label_cdsts",
    "polyA_tail",
    "n"
  )

  stopifnot(
    "df_category must be a data.frame" = is.data.frame(df_category)
  )

  missing_cols <- setdiff(required_cols, colnames(df_category))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "df_category is missing required columns: %s. Check the classification file before calculating proportions.",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (!is.numeric(df_category$n)) {
    stop("df_category$n must be numeric so counts can be summed consistently.")
  }

  if (anyNA(df_category$n)) {
    stop(
      "df_category$n contains NA. Proportion calculations require explicit counts for every row."
    )
  }

  if (any(df_category$n < 0)) {
    stop("df_category$n contains negative values. Counts must be non-negative.")
  }

  for (col_name in c("polyA_tail", "label_premature", "unspliced")) {
    if (!is.logical(df_category[[col_name]])) {
      stop(sprintf("df_category$%s must be logical.", col_name))
    }
    if (anyNA(df_category[[col_name]])) {
      stop(sprintf(
        "df_category$%s contains NA. Resolve missing labels before calculating proportions.",
        col_name
      ))
    }
  }

  for (col_name in c(
    "DHS",
    "sample",
    "strand",
    "label_tts_falling_region",
    "label_premature_termination",
    "label_tts",
    "label_cdsts"
  )) {
    if (anyNA(df_category[[col_name]])) {
      stop(sprintf(
        "df_category$%s contains NA. Resolve missing annotations before calculating proportions.",
        col_name
      ))
    }
  }

  invalid_strand <- unique(df_category$strand[
    !df_category$strand %in% c("+", "-")
  ])
  if (length(invalid_strand) > 0) {
    stop(sprintf(
      "df_category$strand must contain only '+' or '-'. Found: %s",
      paste(invalid_strand, collapse = ", ")
    ))
  }
}


# ---- curate_proportion_categories -------------------------------------------
# Map each read's 3' labels onto one of the feature levels above, with a pA+/- suffix
# recording whether the read carried a polyA tail.
#
# The decision order below is significant: single-exon reads are assigned before the
# premature-termination branch can claim them, and the "before_TTS" cases are tested from
# the most specific label combination to the least, so a read lands in the first matching
# level rather than an arbitrary one.
#
# Parameters:
#   df_category -- data.frame carrying the label columns listed in
#                  validate_proportion_input(); see 00_config.R for the file this table
#                  comes from
#
# Returns:
#   df_category with two added columns:
#     label_3prime     -- the coarse class (single_exon, first_exon, ...)
#     label_3prime_pA  -- class with "_pA+" / "_pA-" appended
curate_proportion_categories <- function(df_category) {
  validate_proportion_input(df_category)

  df_category_curated <- df_category %>%
    mutate(
      IsEnterCDS = case_when(
        strand == "+" ~ mean_TTS >= typical_cdsss,
        strand == "-" ~ mean_TTS <= typical_cdsss,
        TRUE ~ FALSE
      )
    ) %>%
    mutate(
      label_3prime = case_when(
        label_tts_falling_region == "upstream" ~ "upstream",
        label_premature_termination == "singleExon" ~ "single_exon",
        label_tts == "before_TTS" & label_premature ~ "premature_termination",
        label_tts == "before_TTS" &
          label_premature_termination == "firstExon" &
          unspliced ~ "first_exon",
        label_tts == "before_TTS" &
          label_cdsts == "before_CDS" &
          label_premature_termination == "others" ~ "gradual_degradation",
        label_tts == "before_TTS" &
          (label_cdsts == "after_CDS" | label_cdsts == "at_CDS") ~ "UTR3",
        label_tts == "typical_TTS" ~ "typical_TTS",
        label_tts == "overlong_TTS" ~ "over_elongation",
        TRUE ~ "Others"
      )
    ) %>%
    mutate(
      label_3prime_pA = paste0(
        label_3prime,
        "_",
        if_else(polyA_tail, "pA+", "pA-")
      )
    )

  return(df_category_curated)
}


# ---- build_gene_count_long --------------------------------------------------
# Collapse the curated reads into one row per DHS x sample x feature level.
#
# Reads falling upstream of the gene are dropped, since they are not assigned a 3' class.
# Feature levels outside SELECTED_FEATURE_LEVELS are folded into "others", so the level set
# is exactly the denominator contract above.
#
# The ZC3H4 replicates are dropped here: their scale factors are computed and validated,
# but their reads are not part of this table.
#
# Parameters:
#   df_category_curated -- output of curate_proportion_categories()
#
# Returns:
#   data.frame with columns DHS, sample, geno, selected_features, n; one row per
#   DHS x sample x feature level, with absent combinations filled as n = 0.
build_gene_count_long <- function(df_category_curated) {
  df_gene_count_long <- df_category_curated %>%
    filter(label_3prime != "upstream") %>%
    transmute(
      DHS,
      sample,
      geno = extract_genotype_from_sample(sample),
      selected_features = case_when(
        label_3prime_pA %in%
          setdiff(SELECTED_FEATURE_LEVELS, "others") ~ label_3prime_pA,
        TRUE ~ "others"
      ),
      n
    ) %>%
    filter(geno != "ZC3H4") %>%
    mutate(
      selected_features = factor(
        selected_features,
        levels = SELECTED_FEATURE_LEVELS
      ),
      geno = factor(geno, levels = c("WT", "ARS2", "RBM7", "ZC3H18", "ZCCHC8"))
    ) %>%
    group_by(DHS, sample, geno, selected_features) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    tidyr::complete(
      tidyr::nesting(DHS, sample, geno),
      selected_features,
      fill = list(n = 0)
    )

  return(df_gene_count_long)
}


# ---- make_feature_id -------------------------------------------------------
# Convert a readable feature level into its identifier form, for use as a column key:
# "-" becomes "n" and "+" becomes "p".
#
# Parameters:
#   selected_features -- character or factor, e.g. "first_exon_pA+"
#
# Returns:
#   character, e.g. "first_exon_pAp"
make_feature_id <- function(selected_features) {
  feature_id <- selected_features %>%
    as.character() %>%
    str_replace("\\-", "n") %>%
    str_replace("\\+", "p")

  return(feature_id)
}


# ---- calculate_pooled_counts_then_proportion -------------------------------
# Rescale reads by scale factor, pool them per DHS x genotype, and express each feature as
# a proportion of its group's pooled total.
#
# This is the function whose df_freq_long is the heatmap data table.
#
# Parameters:
#   df_category         -- data.frame with columns:
#                            DHS               -- DHS id
#                            sample            -- replicate name
#                            selected_features -- feature level, e.g. "UTR3_pA-"
#                            n                 -- read count for that cell
#                          plus the label columns validate_proportion_input() requires
#   df_scale_factor     -- data.frame with columns `sample` and `scale_factor`, from
#                          define_scale_factor() in 04-Normalization/scale_factors.R
#
# Returns:
#   list with three data.frames:
#     df_category_curated -- the curated per-read table
#     df_gene_count_long  -- columns DHS, sample, geno, selected_features, n, scale_factor,
#                            scaled_n
#     df_freq_long        -- columns DHS, feature_id, geno, feature_count, geno_count,
#                            scaled_freq
calculate_pooled_counts_then_proportion <- function(
  df_category,
  df_scale_factor
) {
  df_category_curated <- curate_proportion_categories(df_category)

  if (is.null(df_scale_factor)) {
    stop(
      "df_scale_factor is required. Derive it with define_scale_factor() in ",
      "04-Normalization/scale_factors.R."
    )
  }
  validate_scale_factor_table(df_scale_factor)

  df_gene_count_long <- build_gene_count_long(df_category_curated) %>%
    dplyr::left_join(
      df_scale_factor,
      by = c("sample" = "sample")
    ) %>%
    mutate(sample = as.character(sample))

  if (anyNA(df_gene_count_long$scale_factor)) {
    missing_scale_samples <- unique(df_gene_count_long$sample[is.na(
      df_gene_count_long$scale_factor
    )])
    stop(sprintf(
      "Missing scale_factor for sample(s): %s. Check define_scale_factor() output before pooling scaled counts.",
      paste(missing_scale_samples, collapse = ", ")
    ))
  }

  if (
    any(!is.finite(df_gene_count_long$scale_factor)) ||
      any(df_gene_count_long$scale_factor <= 0)
  ) {
    stop(
      "df_gene_count_long contains non-finite or non-positive scale_factor values after joining scale factors."
    )
  }

  df_gene_count_long <- df_gene_count_long %>%
    mutate(scaled_n = n / scale_factor) %>%
    mutate(sample = factor(sample, levels = unique(sample)))

  df_total <- df_gene_count_long %>%
    group_by(DHS, geno) %>%
    summarise(geno_count = sum(scaled_n), .groups = "drop")

  df_freq_long <- df_gene_count_long %>%
    mutate(
      feature_id = make_feature_id(selected_features)
    ) %>%
    group_by(DHS, feature_id, geno) %>%
    summarise(
      feature_count = sum(scaled_n),
      .groups = "drop"
    ) %>%
    dplyr::left_join(df_total, by = c("DHS", "geno")) %>%
    mutate(
      scaled_freq = if_else(
        geno_count > 0,
        feature_count / geno_count,
        NA_real_
      )
    )

  return(list(
    df_category_curated = df_category_curated,
    df_gene_count_long = df_gene_count_long,
    df_freq_long = df_freq_long
  ))
}

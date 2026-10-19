# eRNA characterization -- the four relations the manuscript reports
#
# This folder holds the three stages that produce the axes of those panels. The panels themselves are drawn elsewhere; these are the quantities behind them.
#
#   1. eRNA length on the plus strand vs eRNA length on the minus strand
#        both axes from the per-strand representative read-length summary, an input assembled outside this deposit
#   2. log10 STARR-seq signal vs representative eRNA length (both strands)
#        signal from starrseq_per_dhs.R; length from the same read-length summary
#   3. log10 STARR-seq signal vs eRNA TSS spacing
#        signal from starrseq_per_dhs.R; spacing from the TSS-space table written by tf_density.R
#   4. log10 TF peak density vs eRNA TSS spacing
#        density from tf_density.R, which counts ReMap peak blocks with remap_tf_matrix.R; spacing from the same TSS-space table
#
# remap_tf_matrix.R -- count ReMap TF ChIP-seq peak overlaps per DHS region
# Requires GenomicRanges, rtracklayer and stringr; functions only.
# Each path defaults to the file or directory named by an environment variable, so no machine-specific location is fixed here.

library(GenomicRanges)
library(rtracklayer)
library(stringr)

# Shared libraries: dplyr and stringr verbs, and the package attaches this file relies on.
source("../../02-DHS/config_and_libraries.R")

# ── Internal helper ────────────────────────────────────────────────────────────

#' Parse DHS string column to GRanges
#'
#' @description
#' Converts a data.frame with a "DHS" column (format "chrN_start-end") to a GRanges object, preserving all other columns as metadata.
#'
#' @param df_dhs data.frame. Must contain a "DHS" column with strings matching the pattern "chrN_start-end" (e.g. "chr1_12345-56789").
#'
#' @return GRanges with seqnames/start/end derived from DHS and all original columns kept as mcols.
parse_dhs_to_granges <- function(df_dhs) {
  stopifnot(
    "df_dhs must be a data.frame" = is.data.frame(df_dhs),
    "df_dhs must have a 'DHS' column" = "DHS" %in% names(df_dhs)
  )

  dhs_pattern <- "(chr[0-9XYMmt]+)_(\\d+)-(\\d+)"

  valid <- stringr::str_detect(df_dhs$DHS, dhs_pattern)
  if (!all(valid)) {
    n_bad <- sum(!valid)
    stop(sprintf(
      "%d DHS string(s) do not match expected format 'chrN_start-end'. Bad examples: %s",
      n_bad,
      paste(head(df_dhs$DHS[!valid], 3), collapse = ", ")
    ))
  }

  df_dhs$seqnames <- stringr::str_extract(df_dhs$DHS, dhs_pattern, group = 1)
  df_dhs$start <- as.integer(stringr::str_extract(
    df_dhs$DHS,
    dhs_pattern,
    group = 2
  ))
  df_dhs$end <- as.integer(stringr::str_extract(
    df_dhs$DHS,
    dhs_pattern,
    group = 3
  ))

  GenomicRanges::makeGRangesFromDataFrame(df_dhs, keep.extra.columns = TRUE)
}

#' Load ReMap TF category table
#'
#' @param categories_file Character. Path to the TF category table: a TSV with a header row and the columns tf_name, category and subcategory.
#'
#' @return data.frame with columns: tf_name, category, subcategory.
load_tf_categories <- function(categories_file) {
  stopifnot(
    "categories_file must be a single character string" = is.character(
      categories_file
    ) &&
      length(categories_file) == 1
  )
  if (!file.exists(categories_file)) {
    stop(sprintf("TF categories file not found: %s", categories_file))
  }
  utils::read.table(
    categories_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
}

# ── Main functions ─────────────────────────────────────────────────────────────

#' Count ReMap TF peak overlaps per DHS region
#'
#' @description
#' For each TF in tf_names, reads the corresponding BigWig file and counts the number of ReMap peak regions (score = 1 blocks) overlapping each DHS.
#' Returns a data.frame suitable for heatmaps and colSums-style summaries.
#'
#' @param dhs_input GRanges OR data.frame.
#'   - GRanges: used directly. If a "DHS" mcols column exists it is used as
#'     the row label; otherwise labels are derived from coordinates.
#'   - data.frame: must have a "DHS" column with strings like "chr1_12345-56789".
#' @param tf_names Character vector. TF names to query. Each name must match a "<TF>.bw" file in bw_dir; a name with no file is skipped with a warning and produces no column.
#' @param bw_dir Character. Directory holding one "<TF>.bw" BigWig file per TF.
#' @param verbose Logical. Print per-TF progress messages. Default: TRUE.
#'
#' @return data.frame with columns:
#'   - DHS       : DHS identifier string
#'   - <TF_name> : integer column per TF — number of peak blocks overlapping that DHS (0 = absent; ≥1 = present, rarely > 1).
#'
#' @details
#' A TF is considered present in a DHS when at least one contiguous block with score > 0 in its BigWig overlaps any part of the DHS. The cell value is the count of such distinct blocks (usually 0 or 1). For binary presence/absence use df_result[, -1] > 0.
#'
#' Only DHS regions are loaded from each BigWig (rtracklayer `which` argument), keeping memory use proportional to the number of DHS regions, not the genome.
#'
#' @examples
#' \dontrun{
#' df_result <- count_remap_tf_overlap(
#'   dhs_input = df_enhancer_DHS,
#'   tf_names  = c("CTCF", "SOX2", "NANOG", "KLF4", "EP300")
#' )
#'
#' # How many DHS have each TF?
#' colSums(df_result[, -1] > 0)
#' }
count_remap_tf_overlap <- function(
  dhs_input,
  tf_names,
  bw_dir,
  verbose = TRUE
) {
  # === Input Validation ===
  stopifnot(
    "tf_names must be a non-empty character vector" = is.character(tf_names) &&
      length(tf_names) > 0,
    "bw_dir must be a single character string" = is.character(bw_dir) &&
      length(bw_dir) == 1,
    "verbose must be a single logical value" = is.logical(verbose) &&
      length(verbose) == 1
  )

  if (!dir.exists(bw_dir)) {
    stop(sprintf("BigWig directory not found: %s", bw_dir))
  }

  # === 1. Normalise dhs_input to GRanges and extract DHS labels ===
  if (inherits(dhs_input, "GRanges")) {
    gr_dhs <- dhs_input
    if ("DHS" %in% names(GenomicRanges::mcols(gr_dhs))) {
      dhs_ids <- as.character(GenomicRanges::mcols(gr_dhs)$DHS)
    } else {
      dhs_ids <- paste0(
        as.character(GenomicRanges::seqnames(gr_dhs)),
        "_",
        GenomicRanges::start(gr_dhs),
        "-",
        GenomicRanges::end(gr_dhs)
      )
    }
  } else if (is.data.frame(dhs_input)) {
    gr_dhs <- parse_dhs_to_granges(dhs_input)
    dhs_ids <- dhs_input$DHS
  } else {
    stop(
      "dhs_input must be a GRanges object or a data.frame with a 'DHS' column."
    )
  }

  n_dhs <- length(gr_dhs)
  if (n_dhs == 0) {
    stop("dhs_input contains 0 regions.")
  }

  # === 2. Validate TF BigWig availability ===
  available_bw <- tools::file_path_sans_ext(
    list.files(bw_dir, pattern = "\\.bw$", full.names = FALSE)
  )

  tf_missing <- setdiff(tf_names, available_bw)
  if (length(tf_missing) > 0) {
    warning(sprintf(
      "%d TF(s) not found as .bw files in '%s' — skipped (no column produced): %s",
      length(tf_missing),
      bw_dir,
      paste(tf_missing, collapse = ", ")
    ))
  }

  tf_to_process <- intersect(tf_names, available_bw)
  # Keep the order the caller asked for
  tf_to_process <- tf_names[tf_names %in% tf_to_process]

  if (length(tf_to_process) == 0) {
    stop(sprintf(
      "No valid TF BigWig files found. Checked %d TF(s) in: %s",
      length(tf_names),
      bw_dir
    ))
  }

  # === 3. Count overlapping peak blocks per TF ===
  if (verbose) {
    message(sprintf(
      "Counting overlaps for %d TF(s) across %d DHS regions...",
      length(tf_to_process),
      n_dhs
    ))
  }

  result_matrix <- matrix(
    0L,
    nrow = n_dhs,
    ncol = length(tf_to_process),
    dimnames = list(NULL, tf_to_process)
  )

  for (i in seq_along(tf_to_process)) {
    tf <- tf_to_process[i]
    bw_path <- file.path(bw_dir, paste0(tf, ".bw"))

    if (verbose) {
      message(sprintf("  [%d/%d] %s", i, length(tf_to_process), tf))
    }

    # Import only the DHS windows (avoids loading the whole genome into memory)
    tf_peaks <- tryCatch(
      rtracklayer::import.bw(bw_path, which = gr_dhs),
      error = function(e) {
        warning(sprintf(
          "Failed to import BigWig for TF '%s': %s",
          tf,
          conditionMessage(e)
        ))
        NULL
      }
    )

    if (is.null(tf_peaks) || length(tf_peaks) == 0) {
      next
    }

    # The BigWig may contain zero-score padding; keep only actual peaks
    tf_peaks <- tf_peaks[GenomicRanges::score(tf_peaks) > 0]
    if (length(tf_peaks) == 0) {
      next
    }

    # Count distinct peak blocks overlapping each DHS
    hits <- GenomicRanges::findOverlaps(gr_dhs, tf_peaks)
    n_hits_per_dhs <- tabulate(S4Vectors::queryHits(hits), nbins = n_dhs)
    result_matrix[, tf] <- as.integer(n_hits_per_dhs)
  }

  # === 4. Assemble output ===
  df_result <- data.frame(DHS = dhs_ids, result_matrix, check.names = FALSE)
  df_result
}

#' Summarise TF overlap counts by category or subcategory
#'
#' @description
#' Given the output of count_remap_tf_overlap(), groups TF columns by their category (REGULATORY / COFACTOR / NON_REGULATORY) or subcategory (e.g. Zinc_Finger, Polycomb) and sums n_hit values within each group.
#'
#' @param df_tf_overlap data.frame. Output of count_remap_tf_overlap():
#'   first column "DHS", remaining columns are TF n_hit counts.
#' @param categories_file Character. Path to the TF category table: a TSV with a header row and the columns tf_name, category and subcategory.
#' @param by Character. Grouping level: "category" or "subcategory".
#'   Default: "category".
#'
#' @return data.frame with one row per DHS and one column per group, where values are the sum of n_hit across all TFs in that group.
#'
#' @examples
#' \dontrun{
#' df_by_cat    <- summarise_tf_by_category(df_result)
#' df_by_subcat <- summarise_tf_by_category(df_result, by = "subcategory")
#' colSums(df_by_cat[, -1])
#' }
summarise_tf_by_category <- function(
  df_tf_overlap,
  categories_file,
  by = c("category", "subcategory")
) {
  by <- match.arg(by)

  stopifnot(
    "df_tf_overlap must be a data.frame" = is.data.frame(df_tf_overlap),
    "'DHS' column must exist in df_tf_overlap" = "DHS" %in% names(df_tf_overlap)
  )

  tf_cols <- setdiff(names(df_tf_overlap), "DHS")
  if (length(tf_cols) == 0) {
    stop(
      "No TF columns found in df_tf_overlap (expected columns beyond 'DHS')."
    )
  }

  df_cats <- load_tf_categories(categories_file)

  # Warn about TF columns not present in the categories table
  unmatched <- setdiff(tf_cols, df_cats$tf_name)
  if (length(unmatched) > 0) {
    warning(sprintf(
      "%d TF column(s) not found in categories file and will be excluded: %s",
      length(unmatched),
      paste(unmatched, collapse = ", ")
    ))
  }

  df_cats_sub <- df_cats[df_cats$tf_name %in% tf_cols, ]
  group_labels <- df_cats_sub[[by]]
  groups <- unique(group_labels)

  result_mat <- vapply(
    groups,
    function(grp) {
      tfs_in_grp <- df_cats_sub$tf_name[group_labels == grp]
      tfs_in_grp <- intersect(tfs_in_grp, tf_cols)
      if (length(tfs_in_grp) == 0) {
        return(rep(0L, nrow(df_tf_overlap)))
      }
      if (length(tfs_in_grp) == 1) {
        return(as.integer(df_tf_overlap[[tfs_in_grp]]))
      }
      as.integer(rowSums(df_tf_overlap[, tfs_in_grp, drop = FALSE]))
    },
    integer(nrow(df_tf_overlap))
  )

  df_out <- data.frame(DHS = df_tf_overlap$DHS, result_mat, check.names = FALSE)
  names(df_out)[-1] <- groups
  df_out
}

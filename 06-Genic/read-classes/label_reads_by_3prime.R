# label_reads_by_3prime.R -- classify each read's 3' end relative to its gene
# Reads the raw BED12 reads for each DHS x sample with read_TLDR_reads_of_one_DHS_in_one_sample() from 02-DHS/raw_reads_access.R.

# ---- standardize_read_format -----------------------------------------------
# Convert V1-V14 raw BED12 format to standardized readable column names, or validate an already-standardized table (one with a "chr" column).
#
# V1-V14 mapping:
#   V1=chr, V2=start, V3=end, V4=readID, V5=polyA_length, V6=strand, V10=exon_count, V13=DHS, V14=sample
# TSS/TTS are computed from strand direction.
standardize_read_format <- function(df_reads) {
  stopifnot(
    "df_reads must be data.frame" = is.data.frame(df_reads),
    "df_reads cannot be empty" = nrow(df_reads) > 0
  )

  # Already transformed?
  if ("chr" %in% names(df_reads)) {
    required <- c(
      "chr",
      "start",
      "end",
      "readID",
      "strand",
      "DHS",
      "sample",
      "polyA_length",
      "exon_count"
    )
    missing <- setdiff(required, names(df_reads))
    if (length(missing) > 0) {
      stop(sprintf(
        "Transformed format incomplete. Missing: %s",
        paste(missing, collapse = ", ")
      ))
    }
    if (!all(c("TSS", "TTS") %in% names(df_reads))) {
      df_reads <- df_reads %>%
        mutate(
          TSS = ifelse(strand == "+", start, end),
          TTS = ifelse(strand == "+", end, start)
        )
    }
    return(df_reads)
  }

  # V1-V14 format?
  if ("V1" %in% names(df_reads)) {
    required_v <- paste0("V", c(1:6, 10, 13, 14))
    missing_v <- setdiff(required_v, names(df_reads))
    if (length(missing_v) > 0) {
      stop(sprintf(
        "V1-V14 format incomplete. Missing: %s",
        paste(missing_v, collapse = ", ")
      ))
    }
    df_reads <- df_reads %>%
      transmute(
        chr = V1,
        start = V2,
        end = V3,
        readID = V4,
        polyA_length = V5,
        strand = V6,
        score1 = V7,
        score2 = V8,
        metadata = V9,
        exon_count = V10,
        exon_lengths = V11,
        strand2 = V12,
        DHS = V13,
        sample = V14,
        TSS = ifelse(strand == "+", start, end),
        TTS = ifelse(strand == "+", end, start)
      )
    return(df_reads)
  }

  stop(sprintf(
    "Unrecognized read format. Expected 'chr' (transformed) or 'V1' (raw V1-V14) column.\nFound: %s",
    paste(names(df_reads), collapse = ", ")
  ))
}

# ---- classify_position_relative_to_gene_structure --------------------------
# Vectorized classification of positions as exonic, intronic, upstream, or downstream relative to a gene's exon structure.
#
# Returns character vector: "exon_N", "intron_N", "upstream", "downstream", or NA.
classify_position_relative_to_gene_structure <- function(
  positions,
  exon_starts,
  exon_ends,
  strand
) {
  stopifnot(
    "positions must be numeric" = is.numeric(positions),
    "positions cannot be empty" = length(positions) > 0,
    "strand must be '+' or '-'" = strand %in% c("+", "-")
  )
  if (is.na(strand)) {
    stop("strand cannot be NA")
  }

  exon_struct <- parse_exon_structure(
    exon_starts,
    exon_ends,
    strand,
    validate = TRUE
  )
  if (is.null(exon_struct)) {
    stop(
      "Failed to parse exon structure. Check that exon coordinates are valid."
    )
  }

  n_pos <- length(positions)
  result <- character(n_pos)

  for (i in 1:n_pos) {
    pos <- positions[i]
    if (is.na(pos)) {
      result[i] <- NA_character_
      next
    }

    exon_info <- position_in_exon(pos, exon_struct)
    if (exon_info$in_exon) {
      result[i] <- sprintf("exon_%d", exon_info$exon_num)
      next
    }

    if (position_is_intronic(pos, exon_struct)) {
      intron_num <- find_intron_number(pos, exon_struct)
      if (is.na(intron_num)) {
        stop(sprintf(
          "Position %d is intronic but find_intron_number returned NA.",
          pos
        ))
      }
      result[i] <- sprintf("intron_%d", intron_num)
      next
    }

    if (strand == "+") {
      result[i] <- if (pos < exon_struct$gene_start) {
        "upstream"
      } else {
        "downstream"
      }
    } else {
      result[i] <- if (pos > exon_struct$gene_end) "upstream" else "downstream"
    }
  }
  return(result)
}

# ---- add_exon_length -------------------------------------------------------
# Vectorized helper: sum the exon sizes from a BED12 exon_sizes field (comma-separated string of exon block lengths).
add_exon_length <- Vectorize(function(x) {
  as.numeric(str_split(x, ",")[[1]]) %>% sum()
})

# ---- get_index_of_tts_vectorized -------------------------------------------
# Classify each TTS position into exon_N, intron_N, upstream, or downstream relative to a gene's exon structure. Vectorized over multiple TTS positions for a single gene.
#
# Parameters:
#   tts        — numeric vector of TTS positions
#   strand     — "+" or "-"
#   exonStarts — comma-separated exon start positions (genomic order)
#   exonEnds   — comma-separated exon end positions (genomic order)
#
# Returns:
#   character vector: "exon_N", "intron_N", "upstream", "downstream"
get_index_of_tts_vectorized <- function(tts, strand, exonStarts, exonEnds) {
  tmp_exonStarts <- exonStarts %>%
    str_split(",") %>%
    .[[1]] %>%
    .[. != ""] %>%
    as.numeric()
  tmp_exonEnds <- exonEnds %>%
    str_split(",") %>%
    .[[1]] %>%
    .[. != ""] %>%
    as.numeric()

  n_exons <- length(tmp_exonStarts)
  n_tts <- length(tts)
  result <- character(n_tts)

  # Vectorized exon check
  tts_matrix <- matrix(rep(tts, each = n_exons), nrow = n_exons)
  starts_matrix <- matrix(rep(tmp_exonStarts, n_tts), nrow = n_exons)
  ends_matrix <- matrix(rep(tmp_exonEnds, n_tts), nrow = n_exons)

  in_exon <- (tts_matrix >= starts_matrix) & (tts_matrix <= ends_matrix)
  exon_indices <- apply(in_exon, 2, function(x) ifelse(any(x), which(x)[1], 0))

  in_exon_mask <- exon_indices > 0
  if (strand == "+") {
    result[in_exon_mask] <- paste0("exon_", exon_indices[in_exon_mask])
  } else {
    result[in_exon_mask] <- paste0(
      "exon_",
      n_exons - exon_indices[in_exon_mask] + 1
    )
  }

  # Check introns for remaining positions
  remaining <- !in_exon_mask
  if (n_exons > 1 && any(remaining)) {
    for (i in 1:(n_exons - 1)) {
      in_intron <- remaining &
        (tts > tmp_exonEnds[i]) &
        (tts < tmp_exonStarts[i + 1])
      if (strand == "+") {
        result[in_intron] <- paste0("intron_", i)
      } else {
        result[in_intron] <- paste0("intron_", n_exons - i)
      }
      remaining <- remaining & !in_intron
    }
  }

  # Upstream / downstream
  if (strand == "+") {
    upstream_mask <- remaining & (tts < tmp_exonStarts[1])
    downstream_mask <- remaining & (tts > tmp_exonEnds[n_exons])
  } else {
    upstream_mask <- remaining & (tts > tmp_exonEnds[n_exons])
    downstream_mask <- remaining & (tts < tmp_exonStarts[1])
  }
  result[upstream_mask] <- "upstream"
  result[downstream_mask] <- "downstream"
  result[result == ""] <- "unknown"

  return(result)
}

# ---- label_reads_by_3prime_categories --------------------------------------
# Full 3' end classification of reads for a single DHS × sample.
#
# Classifies each read's TTS (the 3' end of the read on its gene's strand) relative to:
#   - typical TTS position (label_tts: typical_TTS / before_TTS / overlong_TTS)
#   - CDS end position (label_cdsts: at_CDS / before_CDS / after_CDS)
#   - exon structure (label_tts_falling_region: exon_N / intron_N / upstream / downstream)
#   - premature termination (label_premature_termination: singleExon / firstExon / premature / others)
#   - newly degraded bins (label_newly_degraded: before_TTS_100bp / 200bp / 400bp)
#   - unspliced premature (label_unspliced_premature_termination)
#
# Parameters:
#   query_DHS                       — character, single DHS id
#   query_sample                    — character, single sample id e.g. "mES_WT_rep1"
#   grange_tts                      — GRanges of typical TTS positions (with gene_name metadata)
#   grange_cds                      — GRanges of CDS end positions (with gene_name metadata)
#   df_DHS_gene_map_with_exon_info  — data.frame, DHS-gene map with exon structure
#   dir_dhs_depth                   — character, directory holding one subdirectory per DHS, each holding one subdirectory per sample with that sample's raw BED12 reads
#   postion_window_size             — integer, full window width in bp, applied as ±postion_window_size / 2 around the typical TTS and the CDS end (default 30; callers pass 100, giving ±50 bp)
#
# Returns:
#   data.frame with one row per read, columns:
#     readID, gene_name, DHS, strand, TSS, TTS, length, exonLength,
#     readID, gene_name, DHS, strand, TSS, TTS, length, exonLength, typical_end, typical_cds_end, typical_exon_count, typical_second_exon_start, polyA_length, label_tts, label_newly_degraded, label_tts_falling_region, label_cdsts, exon_count, label_premature_termination, label_unspliced_premature_termination
#   Zero rows when the DHS has no reads in that sample, or when the DHS's gene has no entry in grange_tts or grange_cds.
label_reads_by_3prime_categories <- function(
  query_DHS,
  query_sample,
  grange_tts,
  grange_cds,
  df_DHS_gene_map_with_exon_info,
  dir_dhs_depth,
  postion_window_size = 30
) {
  # Load reads for this DHS × sample
  tmp_reads <- read_TLDR_reads_of_one_DHS_in_one_sample(
    query_DHS = query_DHS,
    query_replicate = query_sample,
    dir_files = dir_dhs_depth
  )

  if (nrow(tmp_reads) == 0) {
    return(data.frame())
  }

  # Look up gene info for this DHS
  tmp_gene_direction <- df_DHS_gene_map_with_exon_info %>%
    filter(DHS == query_DHS) %>%
    .$strand
  tmp_gene_name <- df_DHS_gene_map_with_exon_info %>%
    filter(DHS == query_DHS) %>%
    .$gene_name

  stopifnot(length(tmp_gene_name) == 1)

  # Check gene exists in TTS and CDS GRanges
  if (!tmp_gene_name %in% grange_tts$gene_name) {
    print(paste0("typical tts not exist: ", query_DHS, " in ", query_sample))
    return(data.frame())
  }
  if (length(tmp_gene_direction) != 1) {
    print(paste0("gene direction unknown: ", query_DHS, " in ", query_sample))
    return(data.frame())
  }
  if (!tmp_gene_name %in% grange_cds$gene_name) {
    print(paste0("typical cds not exist: ", query_DHS, " in ", query_sample))
    return(data.frame())
  }

  # Get typical positions
  tmp_gene_typical_end <- end(grange_tts[
    grange_tts$gene_name == tmp_gene_name,
  ])
  tmp_gene_typical_cds_end <- end(grange_cds[
    grange_cds$gene_name == tmp_gene_name,
  ])

  tmp_typical_exon_count <- df_DHS_gene_map_with_exon_info %>%
    filter(DHS == query_DHS) %>%
    .$exonCount
  tmp_pos_first_exon_end <- df_DHS_gene_map_with_exon_info %>%
    filter(DHS == query_DHS) %>%
    .$first_exon_end
  tmp_pos_second_exon_start <- df_DHS_gene_map_with_exon_info %>%
    filter(DHS == query_DHS) %>%
    .$second_exon_start
  tmp_pos_exon_starts <- df_DHS_gene_map_with_exon_info %>%
    filter(DHS == query_DHS) %>%
    .$exonStarts
  tmp_pos_exon_ends <- df_DHS_gene_map_with_exon_info %>%
    filter(DHS == query_DHS) %>%
    .$exonEnds

  # Classify reads
  tmp_reads %>%
    filter(V6 == tmp_gene_direction) %>%
    transmute(
      readID = V4,
      gene_name = tmp_gene_name,
      DHS = query_DHS,
      strand = V6,
      TSS = ifelse(V6 == "+", V2, V3),
      TTS = ifelse(V6 == "+", V3, V2),
      length = abs(TTS - TSS),
      exonLength = add_exon_length(V11),
      typical_end = tmp_gene_typical_end,
      typical_cds_end = tmp_gene_typical_cds_end,
      typical_exon_count = tmp_typical_exon_count,
      typical_second_exon_start = tmp_pos_second_exon_start,
      polyA_length = V5,

      # TTS relative to typical TTS
      label_tts = case_when(
        V6 == "+" &
          TTS > tmp_gene_typical_end + postion_window_size / 2 ~ "overlong_TTS",
        V6 == "+" &
          TTS <= tmp_gene_typical_end + postion_window_size / 2 &
          TTS >= tmp_gene_typical_end - postion_window_size / 2 ~ "typical_TTS",
        V6 == "+" &
          TTS < tmp_gene_typical_end - postion_window_size / 2 ~ "before_TTS",
        V6 == "-" &
          TTS > tmp_gene_typical_end + postion_window_size / 2 ~ "before_TTS",
        V6 == "-" &
          TTS <= tmp_gene_typical_end + postion_window_size / 2 &
          TTS >= tmp_gene_typical_end - postion_window_size / 2 ~ "typical_TTS",
        V6 == "-" &
          TTS < tmp_gene_typical_end - postion_window_size / 2 ~ "overlong_TTS"
      ),

      # Newly degraded bins
      label_newly_degraded = case_when(
        V6 == "+" &
          label_tts == "before_TTS" &
          TTS > tmp_gene_typical_end - 100 ~ "before_TTS_100bp",
        V6 == "+" &
          label_tts == "before_TTS" &
          TTS > tmp_gene_typical_end - 200 ~ "before_TTS_200bp",
        V6 == "+" &
          label_tts == "before_TTS" &
          TTS > tmp_gene_typical_end - 400 ~ "before_TTS_400bp",
        V6 == "-" &
          label_tts == "before_TTS" &
          TTS < tmp_gene_typical_end + 100 ~ "before_TTS_100bp",
        V6 == "-" &
          label_tts == "before_TTS" &
          TTS < tmp_gene_typical_end + 200 ~ "before_TTS_200bp",
        V6 == "-" &
          label_tts == "before_TTS" &
          TTS < tmp_gene_typical_end + 400 ~ "before_TTS_400bp",
        T ~ "Others"
      ),

      # Exon/intron region
      label_tts_falling_region = get_index_of_tts_vectorized(
        TTS,
        tmp_gene_direction,
        tmp_pos_exon_starts,
        tmp_pos_exon_ends
      ),

      # CDS-relative
      label_cdsts = case_when(
        V6 == "+" &
          TTS >
            tmp_gene_typical_cds_end + postion_window_size / 2 ~ "after_CDS",
        V6 == "+" &
          TTS <= tmp_gene_typical_cds_end + postion_window_size / 2 &
          TTS >= tmp_gene_typical_cds_end - postion_window_size / 2 ~ "at_CDS",
        V6 == "+" &
          TTS <
            tmp_gene_typical_cds_end - postion_window_size / 2 ~ "before_CDS",
        V6 == "-" &
          TTS >
            tmp_gene_typical_cds_end + postion_window_size / 2 ~ "before_CDS",
        V6 == "-" &
          TTS <= tmp_gene_typical_cds_end + postion_window_size / 2 &
          TTS >= tmp_gene_typical_cds_end - postion_window_size / 2 ~ "at_CDS",
        V6 == "-" &
          TTS < tmp_gene_typical_cds_end - postion_window_size / 2 ~ "after_CDS"
      ),

      exon_count = V10,

      # Premature termination
      label_premature_termination = case_when(
        tmp_typical_exon_count == 1 ~ "singleExon",
        V6 == "+" & TTS <= tmp_pos_first_exon_end ~ "firstExon",
        V6 == "+" & TTS < tmp_pos_second_exon_start ~ "premature",
        V6 == "+" & TTS >= tmp_pos_second_exon_start ~ "others",
        V6 == "-" & TTS >= tmp_pos_first_exon_end ~ "firstExon",
        V6 == "-" & TTS > tmp_pos_second_exon_start ~ "premature",
        V6 == "-" & TTS <= tmp_pos_second_exon_start ~ "others",
        T ~ "ERROR"
      ),

      # Unspliced premature
      label_unspliced_premature_termination = case_when(
        label_premature_termination == "singleExon" ~ "singleExon",
        exon_count == 1 &
          label_premature_termination == "premature" ~ "premature_unspliced",
        T ~ "others"
      )
    )
}

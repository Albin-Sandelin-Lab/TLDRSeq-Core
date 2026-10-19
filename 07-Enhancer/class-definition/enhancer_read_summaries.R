# enhancer_read_summaries.R -- summarise read geometry for enhancer DHS regions

# ---- add_exon_length -------------------------------------------------------
# Sum a BED12 read's exon block sizes (column V11, comma-separated).
#
# Parameters:
#   x — character, single comma-separated exon-size string, the trailing-comma form used by BED12 blockSizes.
#
# Returns:
#   numeric, the summed exonic length.
add_exon_length <- Vectorize(function(x) {
  as.numeric(str_split(x, ",")[[1]]) %>%
    sum()
})

# ---- label_enhancer_reads --------------------------------------------------
# Summarise the reads of one enhancer DHS in one sample. Reads are grouped by strand, unspliced status and polyA-tail presence, and summarised by TSS/TTS quantiles, read length, exonic length, exon count and polyA-tail length. `freq` is each group's share of the DHS's reads.
#
# Parameters:
#   query_DHS    — character, single DHS id; matched against the V13 column of the read table and the DHS directory under `dir_files`.
#   query_sample — character, single sample id; matched against the V14 column of the read table and the sample directory under `dir_files`.
#   dir_files    — character, directory holding one subdirectory per DHS, each with one subdirectory per sample containing that DHS x sample's reads as a BED12-style table.
#
# Returns:
#   data.frame with one row per DHS x sample x strand x unspliced x polyA_tail group, or an empty data.frame when the DHS has no reads in that sample.
#   Columns: DHS, sample, strand, unspliced, polyA_tail, median_TSS, median_TTS, quantile10_TSS, quantile90_TSS, quantile10_TTS, quantile90_TTS, mean_TSS, mean_TTS, median_length, mean_length, median_exonLength, mean_exonLength, mean_exon_count, median_exon_count, median_polyA_tail_length, mean_polyA_tail_length, n, freq.
label_enhancer_reads <- function(query_DHS, query_sample, dir_files) {
  tmp_reads <- read_TLDR_reads_of_one_DHS_in_one_sample(
    query_DHS = query_DHS,
    query_replicate = query_sample,
    dir_files = dir_files
  )

  if (nrow(tmp_reads) == 0) {
    return(data.frame())
  }

  tmp_reads <- tmp_reads %>%
    transmute(
      DHS = V13,
      sample = V14,
      strand = V6,
      TSS = ifelse(V6 == "+", V2, V3),
      TTS = ifelse(V6 == "+", V3, V2),
      length = abs(TTS - TSS),
      exonLength = add_exon_length(V11),
      polyA_length = V5,
      exon_count = V10,
      unspliced = exon_count == 1,
      polyA_tail = polyA_length > 15
    )

  tmp_reads %>%
    group_by(DHS, sample, strand, unspliced, polyA_tail) %>%
    summarise(
      median_TSS = median(TSS),
      median_TTS = median(TTS),
      quantile10_TSS = quantile(TSS, 0.1),
      quantile90_TSS = quantile(TSS, 0.9),
      quantile10_TTS = quantile(TTS, 0.1),
      quantile90_TTS = quantile(TTS, 0.9),
      mean_TSS = mean(TSS),
      mean_TTS = mean(TTS),
      median_length = median(length),
      mean_length = mean(length),
      median_exonLength = median(exonLength),
      mean_exonLength = mean(exonLength),
      mean_exon_count = mean(exon_count),
      median_exon_count = median(exon_count),
      median_polyA_tail_length = median(polyA_length),
      mean_polyA_tail_length = mean(polyA_length),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(freq = n / sum(n))
}

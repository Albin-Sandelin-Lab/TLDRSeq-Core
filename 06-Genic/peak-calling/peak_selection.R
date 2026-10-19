# peak_selection.R -- choose one reference TTS per gene from the called peaks

# ---- pick_most_likely_TTS_peak ---------------------------------------------
# Take the peak positions find_peaks() called for one gene and pick the one that becomes the gene's reference 3' end. Peaks below count_threshold, or below prop_threshold of total_counts, are dropped first. Three candidates are then formed among the survivors and the returned label says which of them coincide, in this order of preference: the peak closest to Refseq annotation, the peak with the most counts, and the most 3' extended in the direction of transcription.
#
# Parameters:
#   strand              — character, "+" or "-", the gene's strand; the plus branch keeps the last element of each candidate set and the minus branch the first
#   local_peaks         — character, comma-separated peak positions on the coordinate scale of refseq_tts_position, strictly increasing; NA when the gene has no called peak
#   local_counts        — character, comma-separated read counts, one per entry of local_peaks
#   refseq_tts_position — numeric, the annotated gene 3' end
#   total_counts        — numeric, total counts over the gene's windows, the denominator of prop_threshold
#   count_threshold     — numeric, minimum peak count
#   prop_threshold      — numeric, minimum peak count as a share of total_counts
#
# Returns:
#   data.frame with columns label_peak_picked (which rule chose the peak) and peak_picked (the chosen position).
pick_most_likely_TTS_peak <- function(
  strand = "-",
  local_peaks,
  local_counts,
  refseq_tts_position,
  total_counts,
  count_threshold = 10,
  prop_threshold = 0.1
) {
  if (is.na(local_peaks)) {
    return(data.frame(
      label_peak_picked = "refseq_TTS",
      peak_picked = refseq_tts_position %>% as.numeric()
    ))
  }

  if (!str_detect(local_peaks, ",")) {
    return(data.frame(
      label_peak_picked = "sole peak",
      peak_picked = local_peaks %>% as.numeric()
    ))
  }

  tmp_peaks <- str_split(local_peaks, ",")[[1]] %>% as.numeric()
  tmp_counts <- str_split(local_counts, ",")[[1]] %>% as.numeric()

  stopifnot(rank(tmp_peaks) == seq(length(tmp_peaks)))

  ## pre-filter
  tmp_indices <- seq(length(tmp_peaks))
  tmp_indices_filtered <- tmp_indices[
    tmp_counts[tmp_indices] > count_threshold &
      tmp_counts[tmp_indices] / total_counts > prop_threshold
  ]

  if (length(tmp_indices_filtered) == 0) {
    return(data.frame(
      label_peak_picked = "threshold failed",
      peak_picked = refseq_tts_position %>% as.numeric()
    ))
  }

  tmp_peaks <- tmp_peaks[tmp_indices_filtered]
  tmp_counts <- tmp_counts[tmp_indices_filtered]
  tmp_indices <- seq(length(tmp_peaks))

  tmp_i_closest_to_refseq_position <- tmp_indices[
    abs(tmp_peaks[tmp_indices] - refseq_tts_position) ==
      min(abs(tmp_peaks[tmp_indices] - refseq_tts_position))
  ]

  tmp_i_peak_with_highest_counts <- tmp_indices[
    tmp_counts[tmp_indices] == max(tmp_counts[tmp_indices])
  ]

  if (strand == "+") {
    tmp_i_closest_to_refseq_position <- tmp_i_closest_to_refseq_position %>%
      tail(1)
    tmp_i_peak_with_highest_counts <- tmp_i_peak_with_highest_counts %>% tail(1)
    tmp_i_peak_elongated_farthest <- tmp_indices %>% tail(1)
  } else {
    tmp_i_closest_to_refseq_position <- tmp_i_closest_to_refseq_position %>%
      head(1)
    tmp_i_peak_with_highest_counts <- tmp_i_peak_with_highest_counts %>% head(1)
    tmp_i_peak_elongated_farthest <- tmp_indices %>% head(1)
  }

  if (
    tmp_i_closest_to_refseq_position == tmp_i_peak_with_highest_counts &
      tmp_i_peak_with_highest_counts == tmp_i_peak_elongated_farthest
  ) {
    return(data.frame(
      label_peak_picked = "max_counts, closest_to_refseq, most_elongated",
      peak_picked = tmp_peaks[tmp_i_closest_to_refseq_position]
    ))
  }

  if (tmp_i_closest_to_refseq_position == tmp_i_peak_with_highest_counts) {
    return(data.frame(
      label_peak_picked = "max_counts, closest_to_refseq",
      peak_picked = tmp_peaks[tmp_i_peak_with_highest_counts]
    ))
  }

  if (tmp_i_closest_to_refseq_position == tmp_i_peak_elongated_farthest) {
    return(data.frame(
      label_peak_picked = "closest_to_refseq, most_elongated",
      peak_picked = tmp_peaks[tmp_i_peak_elongated_farthest]
    ))
  }

  if (tmp_i_peak_elongated_farthest == tmp_i_peak_with_highest_counts) {
    return(data.frame(
      label_peak_picked = "max_counts, farthest_peak",
      peak_picked = tmp_peaks[tmp_i_peak_with_highest_counts]
    ))
  }

  return(data.frame(
    label_peak_picked = "max_counts",
    peak_picked = tmp_peaks[tmp_i_peak_with_highest_counts]
  ))
}

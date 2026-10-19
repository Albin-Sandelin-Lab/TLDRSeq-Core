# dhs_gc_content.R -- GC content and width of each DHS

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  # memes supplies get_sequence(), the range-to-sequence reader extractSequencesFromGR() calls.
  library(memes)
  # Supplies BSgenome.Mmusculus.UCSC.mm39, the default genome of extractSequencesFromGR().
  library(BSgenome.Mmusculus.UCSC.mm39)
})

# ---- extractSequencesFromGR ------------------------------------------------
# Sequence of each range, read from a genome with get_sequence().
#
# Each range is widened by 2 * flank bases and recentred before it is read, so flank = 0 reads the range as it stands and a positive flank reads the range plus that much flanking sequence on each side.
#
# Parameters:
#   gr       — GRanges. Ranges to read; their names are kept on the returned sequences.
#   flank    — integer, bp of flanking sequence added on each side of every range (default 0).
#   bsgenome — BSgenome object for the assembly the ranges are in.
#
# Returns:
#   DNAStringSet, one sequence per range, each of width width(gr) + 2 * flank bp, named as gr.
extractSequencesFromGR <- function(
  gr,
  flank = 0,
  bsgenome = BSgenome.Mmusculus.UCSC.mm39
) {
  gr %>%
    resize(width(gr) + 2 * flank, fix = "center") %>%
    get_sequence(bsgenome)
}

# ---- calc_gc_percent -------------------------------------------------------
# GC content of each sequence, as a percentage of its length.
#
# Every base counts towards the length, including any N, so a sequence of pure N scores 0 rather than missing.
#
# Parameters:
#   dna_set — DNAStringSet. Sequences to measure.
#
# Returns:
#   numeric vector, one GC percentage (0-100) per sequence, in the input order.
calc_gc_percent <- function(dna_set) {
  gc_counts <- letterFrequency(dna_set, letters = c("G", "C"))
  rowSums(gc_counts) / width(dna_set) * 100
}

# ---- get_DHS_gc_content ----------------------------------------------------
# GC content and width of each DHS.
#
# Parameters:
#   granges_DHS — GRanges of DHS intervals, one range per DHS, whose names are the DHS ids (for example "chr1_12345-12789"); the DHS column of the result is taken from those names.
#   bsgenome    — BSgenome object for the assembly the DHS intervals are in.
#
# Returns:
#   data.frame with one row per DHS and columns:
#     DHS    — character, the DHS id
#     gc     — numeric, GC percentage of the DHS sequence (0-100)
#     length — numeric, DHS width in bp
get_DHS_gc_content <- function(granges_DHS, bsgenome) {
  df_DHS_sequences <- extractSequencesFromGR(
    gr = granges_DHS,
    flank = 0,
    bsgenome = bsgenome
  )

  data.frame(
    DHS = names(df_DHS_sequences),
    gc = calc_gc_percent(df_DHS_sequences),
    length = width(df_DHS_sequences)
  )
}

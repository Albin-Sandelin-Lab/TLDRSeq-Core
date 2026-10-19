# conservation.R -- mean phastCons conservation score of DHS, exons and introns
#
# Conservation is read from a phastCons track, whose score at a base is the posterior probability that the base is conserved, on a 0-1 scale.
# The exon and intron ranges come from a stored element table, one GRanges list of exons and introns for the RefSeq Select transcripts; the code that builds that table belongs to the annotation stage and is not included here.
# The element id is taken from the element ranges' own exon_number and intron_number metadata, whose values are written as "Exon 1", "Intron 1" and so on; the regression selects element_id "Exon 1" and "Intron 1" from these tables.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicScores)
})

# ---- get_conservation_of_gr ------------------------------------------------
# Mean phastCons score of each range.
#
# The mean is the sum of the scores over the range's bases divided by the number of bases the range spans, so a base the track leaves unscored counts in the denominator and contributes nothing to the sum.
#
# Parameters:
#   gr       — GRanges. Ranges to score.
#   gsco_hmm — GScores object holding the scores, e.g. from gscores() over a phastCons conservation track for the mm39 assembly.
#
# Returns:
#   numeric vector, one mean phastCons score per range, on the 0-1 scale of the track, in the input order.
get_conservation_of_gr <- function(gr, gsco_hmm) {
  sc_hmm <- gscores(gsco_hmm, gr, summaryFun = function(x) {
    sum(x, na.rm = T) / length(x)
  })
  sc_hmm$default
}

# ---- get_DHS_conservation --------------------------------------------------
# Mean phastCons score of each DHS.
#
# Parameters:
#   granges_DHS — GRanges of DHS intervals, one range per DHS, carrying a DHS metadata column with the DHS ids
#   gsco_hmm    — GScores object holding the scores, e.g. from gscores() over a phastCons conservation track for the mm39 assembly
#
# Returns:
#   data.frame with one row per DHS and columns:
#     DHS          — character, the DHS id
#     conservation — numeric, mean phastCons score of the DHS bases (0-1)
get_DHS_conservation <- function(granges_DHS, gsco_hmm) {
  data.frame(
    DHS = granges_DHS$DHS,
    conservation = get_conservation_of_gr(granges_DHS, gsco_hmm)
  )
}

# ---- get_exon_conservation -------------------------------------------------
# Mean phastCons score of each exon.
#
# Parameters:
#   granges_exons — GRanges of exon intervals carrying metadata columns gene_name and exon_number, whose values are written as "Exon 1", "Exon 2" and so on
#   gsco_hmm      — GScores object holding the scores, e.g. from gscores() over a phastCons conservation track for the mm39 assembly
#
# Returns:
#   data.frame with one row per exon and columns:
#     gene_name      — character, gene the exon belongs to
#     element_length — numeric, exon length in bp
#     element_id     — character, the exon number, e.g. "Exon 1"
#     conservation   — numeric, mean phastCons score of the exon bases (0-1)
get_exon_conservation <- function(granges_exons, gsco_hmm) {
  granges_exons$element_id <- granges_exons$exon_number

  data.frame(
    gene_name = granges_exons$gene_name,
    element_length = width(granges_exons),
    element_id = granges_exons$element_id,
    conservation = get_conservation_of_gr(granges_exons, gsco_hmm)
  )
}

# ---- get_intron_conservation -----------------------------------------------
# Mean phastCons score of each intron.
#
# Parameters:
#   granges_introns — GRanges of intron intervals carrying metadata columns gene_name and intron_number, whose values are written as "Intron 1", "Intron 2" and so on
#   gsco_hmm        — GScores object holding the scores, e.g. from gscores() over a phastCons conservation track for the mm39 assembly
#
# Returns:
#   data.frame with one row per intron and columns:
#     gene_name      — character, gene the intron belongs to
#     element_length — numeric, intron length in bp
#     element_id     — character, the intron number, e.g. "Intron 1"
#     conservation   — numeric, mean phastCons score of the intron bases (0-1)
get_intron_conservation <- function(granges_introns, gsco_hmm) {
  granges_introns$element_id <- granges_introns$intron_number

  data.frame(
    gene_name = granges_introns$gene_name,
    element_length = width(granges_introns),
    element_id = granges_introns$element_id,
    conservation = get_conservation_of_gr(granges_introns, gsco_hmm)
  )
}

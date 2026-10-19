# isochore_gc.R -- GC content of the isochores overlapping each DHS
#
# An isochore is a long genomic region of roughly uniform GC content, so a DHS takes the mean GC content of the isochores it overlaps rather than the GC content of its own bases.
# The gene-centric regression reads this table as the GC_isochore predictor.
# The value keeps the unit of the isochore table's avg_gc column, which this code carries through unchanged.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(S4Vectors)
  library(dplyr)
})

# ---- get_DHS_isochore_map --------------------------------------------------
# Pair each DHS with the isochores it overlaps.
#
# Each DHS is expanded by DHS_expand bp before the overlap search, so DHS_expand = 0 pairs a DHS only with isochores it overlaps as it stands.
# A DHS overlapping several isochores appears once per isochore, and duplicated_DHS marks every row of such a DHS.
#
# Parameters:
#   granges_DHS      — GRanges of DHS intervals, one range per DHS, carrying a DHS metadata column with the DHS ids
#   grange_isochore  — GRanges of isochore intervals, carrying metadata columns avg_gc (numeric, the isochore's GC content) and isochore_type (character, the isochore class)
#   DHS_expand       — integer, bp added on each side of every DHS before the overlap search (default 0)
#
# Returns:
#   data.frame with one row per DHS and isochore overlap, columns:
#     the DHS metadata columns, including DHS
#     the isochore metadata columns, including avg_gc and isochore_type
#     strand         — character, strand of the isochore interval
#     duplicated_DHS — logical, TRUE on every row of a DHS that overlaps more than one isochore
get_DHS_isochore_map <- function(granges_DHS, grange_isochore, DHS_expand = 0) {
  granges_DHS_expand <- granges_DHS + DHS_expand

  tmp_overlaps_DHS_genes <- findOverlapPairs(
    granges_DHS_expand,
    grange_isochore
  )

  tmp_DHS_genes_mapped <- S4Vectors::first(tmp_overlaps_DHS_genes)
  names(tmp_DHS_genes_mapped) <- NULL

  tmp_gene_info <- S4Vectors::second(tmp_overlaps_DHS_genes)
  names(tmp_gene_info) <- NULL

  df_DHS_genes_mapped <- cbind(
    mcols(tmp_DHS_genes_mapped),
    mcols(tmp_gene_info),
    strand = BiocGenerics::strand(tmp_gene_info)
  ) %>%
    as.data.frame() %>%
    mutate(
      duplicated_DHS = (duplicated(DHS) | duplicated(DHS, fromLast = TRUE))
    ) # label those duplicated overlapps

  return(df_DHS_genes_mapped)
}

# ---- summarise_DHS_isochore_gc ---------------------------------------------
# Reduce the DHS and isochore pairs to one mean isochore GC content per DHS.
#
# Parameters:
#   df_DHS_mapped_with_isochores — data.frame. The output of get_DHS_isochore_map().
#   dhs_keep                     — character vector. DHS ids to keep, i.e. the DHS that carry a gene assignment in the DHS-gene map the regression is built on; DHS outside this vector are dropped.
#
# Returns:
#   data.frame with one row per kept DHS and columns:
#     DHS              — character, the DHS id
#     avg_gc           — numeric, mean of the avg_gc values of the isochores overlapping the DHS, in the unit of the isochore table
#     n                — integer, number of isochores the mean was taken over
#     isochore_types   — character, the overlapping isochore types pasted together with ","
summarise_DHS_isochore_gc <- function(df_DHS_mapped_with_isochores, dhs_keep) {
  df_DHS_mapped_with_isochores %>%
    group_by(DHS) %>%
    summarise(
      avg_gc = mean(avg_gc),
      n = n(),
      isochore_types = paste0(isochore_type, collapse = ",")
    ) %>%
    filter(DHS %in% dhs_keep)
}

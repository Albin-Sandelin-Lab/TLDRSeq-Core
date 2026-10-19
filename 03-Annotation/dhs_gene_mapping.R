# dhs_gene_mapping.R -- map DHS regions to genes by TSS overlap

# ---- get_DHS_gene_map ------------------------------------------------------
# Compute DHS-to-gene overlaps from a GRanges of gene TSSs.
# This is the core overlap function used by all annotation generators.
#
# Parameters:
#   grange_gene_tss  — GRanges of gene TSS positions in memory (metadata: gene_id)
#   file_dhs_granges — character, path to an RDS holding a GRanges of DHS regions (metadata: DHS)
#   DHS_expand       — integer, bp to expand each DHS before overlapping (default 0)
#
# Returns:
#   data.frame with columns: DHS, gene_id, ..., duplicated_DHS, duplicated_gene
get_DHS_gene_map <- function(
  grange_gene_tss,
  file_dhs_granges,
  DHS_expand = 0
) {
  granges_DHS <- readRDS(file_dhs_granges)
  granges_DHS_expand <- granges_DHS + DHS_expand

  tmp_overlaps_DHS_genes <- findOverlapPairs(
    granges_DHS_expand,
    grange_gene_tss
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
      duplicated_DHS = (duplicated(DHS) | duplicated(DHS, fromLast = TRUE)),
      duplicated_gene = (duplicated(gene_id) |
        duplicated(gene_id, fromLast = TRUE))
    )

  return(df_DHS_genes_mapped)
}

# codon_usage.R -- codon adaptation index per gene
#
# The codon adaptation index summarises how a coding sequence uses the genetic code, computed from a codon usage table of the coding sequences so that each gene is compared with a reference set of genes.
# The reference set here is the ribosomal-protein genes, whose codon usage is the one the index is measured against; the column is named ribosomal.
# The gene-centric regression reads this column as the CAI predictor.
#
# The coding sequence of a gene is the part of its spliced RNA after the 5' UTR, so the caller slices the sequence set before calling these functions.

suppressPackageStartupMessages({
  library(Biostrings)
  library(coRdon)
  library(dplyr)
  library(stringr)
})

# ---- select_ribosomal_genes ------------------------------------------------
# RefSeq accessions of the ribosomal-protein genes.
#
# These are the genes whose symbol carries "Rp", and they form the reference codon usage set the ribosomal column of calculate_cai() is computed against.
#
# Parameters:
#   df_gene_map — data.frame. Gene table with columns gene_name and refseq_name, e.g. the DHS-gene map with exon structure; one row per DHS, so an accession appears once per row of its gene.
#
# Returns:
#   character vector, the refseq_name of every row whose gene_name carries "Rp".
select_ribosomal_genes <- function(df_gene_map) {
  df_gene_map %>%
    filter(str_detect(gene_name, "Rp")) %>%
    .$refseq_name
}

# ---- calculate_cai ---------------------------------------------------------
# Codon adaptation index of each coding sequence.
#
# Parameters:
#   cds_sequences           — DNAStringSet. Coding sequences to score, named by RefSeq accession; those names become refseq_name.
#   ribosomal_refseq_names  — character vector. Accessions to use as the reference set, e.g. select_ribosomal_genes(df_gene_map); they name members of cds_sequences.
#
# Returns:
#   data.frame with one row per coding sequence and columns:
#     ribosomal   — numeric, the index against the ribosomal-protein codon usage; dimensionless, rising as a sequence's codon usage approaches that of the reference set
#     self        — numeric, the index against the codon usage of the whole input
#     refseq_name — character, the accession
calculate_cai <- function(cds_sequences, ribosomal_refseq_names) {
  tmp_codonTable <- codonTable(cds_sequences)

  tmp_ribo_codonTable <- codonTable(cds_sequences[ribosomal_refseq_names])

  df_cai <- CAI(
    tmp_codonTable,
    subsets = list(ribosomal = tmp_ribo_codonTable, self = tmp_codonTable)
  ) %>%
    as.data.frame() %>%
    mutate(refseq_name = names(cds_sequences))

  return(df_cai)
}

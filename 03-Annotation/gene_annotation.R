# gene_annotation.R -- build gene TSS/TTS/exon annotation from the RefSeq Select tables
# Every builder takes the path of its input annotation table as a required argument; no annotation table is distributed with this repository.

# ---- get_GRange_tss_of_mm39_genes ------------------------------------------
# Build a GRanges of gene TSS positions from a first-3-exon annotation table.
# Each gene's TSS is expanded by window_expand_tss bp on both sides.
#
# Parameters:
#   file_refseq       — character, path to the RDS first-3-exon annotation table (columns: chrom, TSS, strand, gene_name, entrezgene_id, ensembl_gene_id, refseq_name)
#   window_expand_tss — integer, bp to expand around each TSS (default 300)
#
# Returns:
#   GRanges with metadata: gene_id, mgi_symbol, ensembl_gene_id, refseq_name
get_GRange_tss_of_mm39_genes <- function(
  file_refseq,
  window_expand_tss = 300
) {
  granges_tss <- readRDS(file_refseq) %>%
    transmute(
      seqnames = chrom,
      start = TSS,
      end = TSS,
      strand = strand,
      gene_id = entrezgene_id,
      mgi_symbol = gene_name,
      ensembl_gene_id = ensembl_gene_id,
      refseq_name = refseq_name
    ) %>%
    GRanges()

  granges_tss <- granges_tss + window_expand_tss
  return(granges_tss)
}

# ---- get_GRange_tts_of_mm39_genes ------------------------------------------
# Build a GRanges of gene TTS positions from a first-3-exon annotation table.
# Each gene's TTS is expanded by window_expand_tts bp on both sides.
#
# Parameters:
#   file_refseq       — character, path to the RDS first-3-exon annotation table (columns: chrom, TTS, strand, gene_name, entrezgene_id, ensembl_gene_id, refseq_name)
#   window_expand_tts — integer, bp to expand around each TTS (default 300)
#
# Returns:
#   GRanges with metadata: gene_id, mgi_symbol, ensembl_gene_id, refseq_name
get_GRange_tts_of_mm39_genes <- function(
  file_refseq,
  window_expand_tts = 300
) {
  granges_tts <- readRDS(file_refseq) %>%
    transmute(
      seqnames = chrom,
      start = TTS,
      end = TTS,
      strand = strand,
      gene_id = entrezgene_id,
      mgi_symbol = gene_name,
      ensembl_gene_id = ensembl_gene_id,
      refseq_name = refseq_name
    ) %>%
    GRanges()

  granges_tts <- granges_tts + window_expand_tts
  return(granges_tts)
}

# ---- get_GRange_of_mm39_genes_general --------------------------------------
# General-purpose GRanges builder: can produce TSS, TTS, or full gene-body ranges from a first-3-exon annotation table.
#
# Parameters:
#   file_refseq        — character, path to the RDS first-3-exon annotation table
#   output_mode        — character, "TSS", "TTS", or "full"
#   postprocess_expand — integer, bp to expand each range (default 300)
#
# Returns:
#   GRanges with metadata: gene_id, entrez_gene_id, mgi_symbol, ensembl_gene_id, refseq_name
get_GRange_of_mm39_genes_general <- function(
  file_refseq,
  output_mode = "TSS",
  postprocess_expand = 300
) {
  stopifnot(output_mode %in% c("TSS", "TTS", "full"))

  if (output_mode == "TSS") {
    granges_output <- readRDS(file_refseq) %>%
      transmute(
        seqnames = chrom,
        start = TSS,
        end = TSS,
        strand = strand,
        gene_id = gene_name,
        entrez_gene_id = entrezgene_id,
        mgi_symbol = gene_name,
        ensembl_gene_id = ensembl_gene_id,
        refseq_name = refseq_name
      ) %>%
      GRanges()
  } else if (output_mode == "TTS") {
    granges_output <- readRDS(file_refseq) %>%
      transmute(
        seqnames = chrom,
        start = TTS,
        end = TTS,
        strand = strand,
        gene_id = gene_name,
        entrez_gene_id = entrezgene_id,
        mgi_symbol = gene_name,
        ensembl_gene_id = ensembl_gene_id,
        refseq_name = refseq_name
      ) %>%
      GRanges()
  } else if (output_mode == "full") {
    granges_output <- readRDS(file_refseq) %>%
      transmute(
        seqnames = chrom,
        start = txStart,
        end = txEnd,
        strand = strand,
        gene_id = gene_name,
        entrez_gene_id = entrezgene_id,
        mgi_symbol = gene_name,
        ensembl_gene_id = ensembl_gene_id,
        refseq_name = refseq_name
      ) %>%
      GRanges()
  }

  granges_output <- granges_output + postprocess_expand
  return(granges_output)
}

# ---- get_gene_id_convert_from_biomart --------------------------------------
# Query Ensembl BioMart to convert MGI symbols to Ensembl/Entrez IDs.
# Requires internet access and the biomaRt package.
#
# Parameters:
#   query_values — character vector of MGI gene symbols
#   filter_key   — biomaRt filter name (default "mgi_symbol")
#
# Returns:
#   data.frame: ensembl_gene_id, entrezgene_id, mgi_symbol, refseq_mrna
get_gene_id_convert_from_biomart <- function(
  query_values,
  filter_key = "mgi_symbol"
) {
  library(biomaRt)
  mart <- useMart(
    "ENSEMBL_MART_ENSEMBL",
    dataset = "mmusculus_gene_ensembl",
    host = "https://www.ensembl.org"
  )
  getBM(
    attributes = c(
      "ensembl_gene_id",
      "entrezgene_id",
      "mgi_symbol",
      "refseq_mrna"
    ),
    filters = filter_key,
    values = query_values,
    mart = mart
  )
}

# ---- get_gene_annotation_from_refseq_select --------------------------------
# Build a first-3-exon annotation table from a RefSeq Select CSV export.
# Computes TSS, TTS and first/second/third exon boundaries for each transcript, then adds Ensembl/Entrez IDs via BioMart.
#
# Parameters:
#   file_refseq_csv — character, path to the RefSeq Select CSV (BED12-like columns: chrom, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, strand, name, name2)
#
# Returns:
#   data.frame with columns: gene_name, refseq_name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, TSS, TTS, first_exon_start, first_exon_end, second_exon_start, second_exon_end, third_exon_start, third_exon_end, isProteinCoding, ensembl_gene_id, entrezgene_id
get_gene_annotation_from_refseq_select <- function(file_refseq_csv) {
  df_refseq <- read.csv(file_refseq_csv)

  df_first_second_exon <- df_refseq %>%
    mutate(
      TSS = case_when(strand == "+" ~ txStart, strand == "-" ~ txEnd),
      TTS = case_when(strand == "-" ~ txStart, strand == "+" ~ txEnd),
      first_exon_start = case_when(
        exonCount == 1 & strand == "+" ~ str_extract(exonStarts, "(\\d+),", 1),
        exonCount == 1 & strand == "-" ~ str_extract(exonEnds, "(\\d+),", 1),
        strand == "+" ~ str_extract(exonStarts, "(\\d+),", 1),
        strand == "-" ~ str_extract(exonEnds, "[0-9,]*,(\\d+),$", 1),
        T ~ "ERROR"
      ),
      second_exon_start = case_when(
        exonCount == 1 ~ NA,
        exonCount == 2 & strand == "+" ~ str_extract(
          exonStarts,
          "\\d+,(\\d+),",
          1
        ),
        exonCount == 2 & strand == "-" ~ str_extract(
          exonEnds,
          "(\\d+),\\d+,$",
          1
        ),
        strand == "+" ~ str_extract(exonStarts, "\\d+,(\\d+),", 1),
        strand == "-" ~ str_extract(exonEnds, "(\\d+),\\d+,$", 1),
        T ~ "ERROR"
      ),
      third_exon_start = case_when(
        exonCount <= 2 ~ NA,
        strand == "+" ~ str_extract(exonStarts, "\\d+,\\d+,(\\d+),", 1),
        strand == "-" ~ str_extract(exonEnds, "(\\d+),\\d+,\\d+,$", 1),
        T ~ "ERROR"
      ),
      first_exon_end = case_when(
        exonCount == 1 & strand == "+" ~ str_extract(exonEnds, "(\\d+),", 1),
        exonCount == 1 & strand == "-" ~ str_extract(exonStarts, "(\\d+),", 1),
        strand == "+" ~ str_extract(exonEnds, "(\\d+),", 1),
        strand == "-" ~ str_extract(exonStarts, "[0-9,]*,(\\d+),$", 1),
        T ~ "ERROR"
      ),
      second_exon_end = case_when(
        exonCount == 1 ~ NA,
        exonCount == 2 & strand == "+" ~ str_extract(
          exonEnds,
          "\\d+,(\\d+),",
          1
        ),
        exonCount == 2 & strand == "-" ~ str_extract(
          exonStarts,
          "(\\d+),\\d+,$",
          1
        ),
        strand == "+" ~ str_extract(exonEnds, "\\d+,(\\d+),", 1),
        strand == "-" ~ str_extract(exonStarts, "(\\d+),\\d+,$", 1),
        T ~ "ERROR"
      ),
      third_exon_end = case_when(
        exonCount <= 2 ~ NA,
        strand == "+" ~ str_extract(exonEnds, "\\d+,\\d+,(\\d+),", 1),
        strand == "-" ~ str_extract(exonStarts, "(\\d+),\\d+,\\d+,$", 1),
        T ~ "ERROR"
      )
    ) %>%
    mutate(
      first_exon_start = as.numeric(first_exon_start),
      second_exon_start = as.numeric(second_exon_start),
      third_exon_start = as.numeric(third_exon_start),
      first_exon_end = as.numeric(first_exon_end),
      second_exon_end = as.numeric(second_exon_end),
      third_exon_end = as.numeric(third_exon_end)
    ) %>%
    transmute(
      gene_name = name2,
      refseq_name = name,
      chrom,
      strand,
      txStart,
      txEnd,
      cdsStart,
      cdsEnd,
      exonCount,
      exonStarts,
      exonEnds,
      TSS,
      TTS,
      first_exon_start,
      first_exon_end,
      second_exon_start,
      second_exon_end,
      third_exon_start,
      third_exon_end,
      isProteinCoding = str_detect(name, "NM")
    )

  tmp_symbol_map <- get_gene_id_convert_from_biomart(
    query_values = unique(df_first_second_exon$gene_name),
    filter_key = "mgi_symbol"
  ) %>%
    dplyr::select(c("ensembl_gene_id", "entrezgene_id", "mgi_symbol")) %>%
    filter(!duplicated(mgi_symbol))

  df_first_second_exon_annotated <- df_first_second_exon %>%
    dplyr::left_join(tmp_symbol_map, by = c("gene_name" = "mgi_symbol"))

  return(df_first_second_exon_annotated)
}

# ---- get_non_overlap_exonic_region_based_on_refseq_select ------------------
# Build a GRanges of merged non-overlapping exonic regions from a RefSeq Select GTF. Requires GenomicFeatures; the transcript database is built from the supplied GTF, which must be the RefSeq Select track for the same genome assembly as the DHS ranges.
#
# Parameters:
#   file_gtf — character, path to the RefSeq Select GTF
#
# Returns:
#   GRanges of reduced (merged) exonic regions, strand-agnostic
get_non_overlap_exonic_region_based_on_refseq_select <- function(file_gtf) {
  library(GenomicFeatures)
  txdb <- makeTxDbFromGFF(file_gtf)
  GenomicRanges::reduce(exons(txdb), ignore.strand = TRUE)
}

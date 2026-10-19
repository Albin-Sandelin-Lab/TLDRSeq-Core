# exon_intron_gc.R -- GC content of a gene's exons and its introns
#
# Each exon and each intron is measured on its own sequence, so the table carries one GC percentage per element rather than one per gene.
# Elements are numbered in transcript order from the 5' end, which is what gives element "exon 1" and element "intron 1" their meaning.
# The gene-centric regression takes the gc of elements "exon 1" and "intron 1" from this table and derives the Exon1_AT and Intron1_AT predictors from them as 100 minus that gc.
#
# Requires calc_gc_percent() from 06-Genic/covariates/dhs_gc_content.R.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(GenomicFeatures)
  library(dplyr)
  library(stringr)
  library(BSgenome.Mmusculus.UCSC.mm39)
})

# ---- extract_exon_intron_seqs ----------------------------------------------
# The exons and introns of one RefSeq transcript, each with the GC percentage of its own sequence.
#
# Parameters:
#   nm_id  — character scalar. RefSeq mRNA accession with version, which is the transcript's GENEID in the TxDb, e.g. "NM_000000.1".
#   txdb   — TxDb object built from the RefSeq Select GTF for the mm39 assembly.
#   genome — BSgenome object for the mm39 assembly.
#
# Returns:
#   list with:
#     exons      — GRanges of the transcript's exons, with metadata refseq_id and gc (gc = GC percentage of the exon sequence, 0-100)
#     introns    — GRanges of the transcript's introns with the same metadata, or NA when the transcript is single-exon
#     exon_seq   — DNAStringSet of the exon sequences
#     intron_seq — DNAStringSet of the intron sequences, or NA when the transcript is single-exon
#
# The call stops when the accession has no version, when the TxDb holds no transcript for it, or when the transcript carries no transcript id.
extract_exon_intron_seqs <- function(nm_id, txdb, genome) {
  stopifnot(str_detect(nm_id, "NM_\\d+\\.\\d+"))

  # Match against GENEID in txdb (these have versions, so partial match is fine)
  all_geneids <- keys(txdb, keytype = "GENEID")
  geneid_match <- all_geneids[ all_geneids == nm_id ]
  print(geneid_match)
  if (length(geneid_match) == 0) stop("No matching transcript in txdb for: ", nm_id)

  # Get transcript(s) for this NM accession
  txs <- transcriptsBy(txdb, by="gene")[[geneid_match]]
  if (is.null(txs)) stop("No transcripts found for geneid: ", geneid_match)

  tx_name <- mcols(txs)$tx_id
  stopifnot(length(tx_name) != 0)
  # Extract exons & introns
  exons <- exonsBy(txdb, by="tx", use.names=TRUE)[tx_name][[1]]
  introns <- intronsByTranscript(txdb, use.names=TRUE)[tx_name][[1]]

  ## Exon

  exon_seq = getSeq(genome, exons)
  exons$refseq_id = geneid_match
  exons$gc = calc_gc_percent(exon_seq)

  if(length(introns) != 0){
    intron_seq = getSeq(genome, introns)
    introns$refseq_id = geneid_match
    introns$gc = calc_gc_percent(intron_seq)
  }else{
    introns = NA
    intron_seq = NA
  }
  list(
    exons = exons,
    introns = introns,
    exon_seq   = exon_seq,
    intron_seq = intron_seq
  )
}

# ---- get_intron_and_exon_of_one_gene ---------------------------------------
# The exons and introns of one transcript, or NA when they cannot be built.
#
# A transcript the TxDb cannot resolve is reported as a message naming the accession, and NA is returned for it rather than an error propagating.
#
# Parameters:
#   query_gene      — character scalar. RefSeq mRNA accession with version.
#   ref_txdb        — TxDb object built from the RefSeq Select GTF for the mm39 assembly.
#   ref_genome      — BSgenome object for the mm39 assembly.
#
# Returns:
#   the list returned by extract_exon_intron_seqs(), or NA when that call fails.
get_intron_and_exon_of_one_gene <- function(query_gene,
                                            ref_txdb,
                                            ref_genome = BSgenome.Mmusculus.UCSC.mm39) {

  tryCatch({
    extract_exon_intron_seqs(query_gene,
                             txdb = ref_txdb,
                             genome = ref_genome)
  }, error = function(e) {
    message("Error with DHS ", query_gene, ": ", e$message)
    return(NA)
  })
}

# ---- get_intron_and_exon_gc_of_multiple_genes ------------------------------
# GC percentage of the exons and introns of each of several transcripts, as one table.
#
# Elements are put in transcript order from the 5' end, and within a transcript each type is numbered from the 5' end, so element "exon 1" is the 5'-most exon and element "intron 1" the 5'-most intron.
#
# Parameters:
#   query_genes            — character vector. RefSeq mRNA accessions with version.
#   file_refseq_select_gtf — character, path to the RefSeq Select GTF for the mm39 assembly.
#
# Returns:
#   data.frame with one row per exon or intron and columns:
#     DHS                — character, the accession whose element this row is (the query id, kept under the name the analysis gave it)
#     seqnames           — character, chromosome
#     start              — numeric, element start, bp
#     end                — numeric, element end, bp
#     width              — numeric, element length, bp
#     strand             — character, "+" or "-"
#     refseq_id          — character, the accession
#     gc                 — numeric, GC percentage of the element sequence (0-100)
#     type               — factor, "exon" or "intron"
#     order_helper       — numeric, sort key ordering the elements from the 5' end: the start on the plus strand, minus the end on the minus strand
#     x_start            — numeric, bases of the transcript before the element
#     x_end              — numeric, bases of the transcript up to the end of the element
#     index              — integer, position of the element along the transcript
#     type_specific_index — integer, position of the element among the elements of its own type in its transcript
#     element            — character, type and type_specific_index written together, e.g. "exon 1"
#
# A transcript whose exons and introns cannot be built comes back as NA from get_intron_and_exon_of_one_gene() and is not dropped, so reading that entry stops the call instead of skipping the transcript.
# A single-exon transcript contributes its exons only.
get_intron_and_exon_gc_of_multiple_genes <- function(query_genes, file_refseq_select_gtf){
  txdb <- makeTxDbFromGFF(file_refseq_select_gtf, format="gtf")

  tmp_list <- progressr::with_progress({
    p <- progressr::progressor(along = query_genes)

    lapply(query_genes, function(gene) {
      res <- get_intron_and_exon_of_one_gene(gene, txdb, BSgenome.Mmusculus.UCSC.mm39)
      p()  # update progress
      res
    })
  })

  names(tmp_list) = query_genes

  tmp_list_filtered = tmp_list[sapply(tmp_list, function(x){
    if(length(x)>= 1){return(TRUE)}
    if(length(x)==0){return(FALSE)}
    if(is.na(x)){return(FALSE)}
    return(FALSE)
  })]

  if(length(tmp_list_filtered) == 0) return(data.frame())

  df_exon_gc = lapply(names(tmp_list_filtered),
                      function(x){
                        if(length(tmp_list_filtered[[x]]$exons) == 1 && is.na(tmp_list_filtered[[x]]$exons)){
                          return(data.frame())
                        }
                        tmp_list_filtered[[x]]$exons %>%
                          as.data.frame() %>%
                          mutate(DHS = x)
                      }) %>%
    bind_rows() %>%
    transmute(DHS,
              seqnames,
              start,
              end,
              width,
              strand,
              refseq_id,
              gc,
              type = "exon")

  df_intron_gc = lapply(names(tmp_list_filtered),
                        function(x){
                          if(length(tmp_list_filtered[[x]]$introns) == 1 && is.na(tmp_list_filtered[[x]]$introns)){
                            return(data.frame())
                          }
                          tmp_list_filtered[[x]]$introns %>%
                            as.data.frame() %>%
                            mutate(DHS = x)
                        }) %>%
    bind_rows() %>%
    transmute(
      DHS,seqnames,
      start,
      end,
      width,
      strand,
      refseq_id,
      gc,
      type = "intron")

  df_exon_and_intron = bind_rows(df_exon_gc, df_intron_gc) %>%
    mutate(type = factor(type, levels = c("intron", "exon"))) %>%
    mutate(order_helper = ifelse(strand == "+", start, - end)) %>%
    arrange(refseq_id, order_helper) %>%
    group_by(refseq_id) %>%
    mutate(x_start = cumsum(width) - width,
           x_end = cumsum(width),
           index = row_number()) %>%
    ungroup() %>%
    group_by(refseq_id, type) %>%
    mutate(type_specific_index = row_number(),
           element = paste0(type, " ", type_specific_index)
           ) %>%
    ungroup()

  return(df_exon_and_intron)
}

# dhs_enhancer_mapping.R -- define the DHS enhancer class by ordered overlap rules

# The two summarisers below overlap annotation ranges held in memory with the DHS regions read from `file_dhs_granges`. The annotation side arrives already built because the callers subset and re-label it before overlapping; the overlap itself is get_DHS_gene_map() from 03-Annotation/dhs_gene_mapping.R.

# ---- summarise_overlap_DHS_with_X_ignore_direction --------------------------
# Summarise a DHS x annotation overlap ignoring strand: one row per DHS-gene pair with the overlapped transcripts, plus flags for DHSs or annotations that take part in more than one overlap. `strand` is set to "*".
#
# Parameters:
#   grange_x         — GRanges of annotation ranges, with gene_id and refseq_name metadata.
#   label            — character, the annotation-type label carried into `annot_type`.
#   DHS_expand       — integer, bp to expand each DHS before overlapping (default 0).
#   file_dhs_granges — character, path to an RDS holding a GRanges of DHS regions (metadata column: DHS).
#
# Returns:
#   data.frame: DHS, gene_id, n_overlapped_transcript, overlapped_transcript, DHS_map_multi_annot, annot_map_multi_DHS, annot_type, strand.
summarise_overlap_DHS_with_X_ignore_direction = function(grange_x, label = "", DHS_expand = 0, file_dhs_granges){
  tmp_DHS_gene_map = get_DHS_gene_map(grange_gene_tss = grange_x,
                                                   file_dhs_granges = file_dhs_granges,
                                                   DHS_expand = DHS_expand) %>%
    group_by(DHS, gene_id) %>%
    summarise(n_overlapped_transcript = n(),
              overlapped_transcript = paste0(refseq_name, collapse = ","),
              .groups = "drop") %>%
    mutate(DHS_map_multi_annot = (duplicated(DHS) | duplicated(DHS, fromLast = TRUE)),
           annot_map_multi_DHS = (duplicated(gene_id) | duplicated(gene_id, fromLast = TRUE)),
           annot_type = label,
           strand = "*")

  return(tmp_DHS_gene_map)
}


# ---- summarise_overlap_DHS_with_X_directional -------------------------------
# The stranded counterpart of the summariser above: one row per DHS-gene-strand triple, with the duplicate flags computed within each strand.
#
# Parameters:
#   grange_x         — GRanges of annotation ranges, with gene_id and refseq_name metadata.
#   label            — character, the annotation-type label carried into `annot_type`.
#   DHS_expand       — integer, bp to expand each DHS before overlapping (default 0).
#   file_dhs_granges — character, path to an RDS holding a GRanges of DHS regions (metadata column: DHS).
#
# Returns:
#   data.frame: DHS, gene_id, strand, n_overlapped_transcript, overlapped_transcript, DHS_map_multi_annot, annot_map_multi_DHS, annot_type.
summarise_overlap_DHS_with_X_directional = function(grange_x, label = "", DHS_expand = 0, file_dhs_granges){
  tmp_DHS_gene_map = get_DHS_gene_map(grange_gene_tss = grange_x,
                                                   file_dhs_granges = file_dhs_granges,
                                                   DHS_expand = DHS_expand) %>%
    group_by(DHS, gene_id, strand) %>%
    summarise(n_overlapped_transcript = n(),
              overlapped_transcript = paste0(refseq_name, collapse = ","),
              .groups = "drop") %>%
    group_by(strand) %>%
    mutate(DHS_map_multi_annot = (duplicated(DHS) | duplicated(DHS, fromLast = TRUE)),
           annot_map_multi_DHS = (duplicated(gene_id) | duplicated(gene_id, fromLast = TRUE)),
           annot_type = label) %>%
    ungroup()

  return(tmp_DHS_gene_map)
}


# ---- get_DHS_map_to_stuff ---------------------------------------------------
# Build one overlap table per annotation source and return them as a list. The sources are the RefSeq Select gene annotation as TSS / TTS / full body, MGI non-coding RNA, MGI enhancers, and the merged exonic regions of the RefSeq Select annotation. TSS/TTS ranges are single-base points; full-body ranges carry the transcript extent. Every overlap is read through the two summarisers above.
#
# Parameters:
#   file_refseq_select_first_3_exons  — character, path to an RDS data.frame of the first-three-exon RefSeq Select annotation (columns: chrom, TSS, TTS, txStart, txEnd, strand, gene_name, entrezgene_id, ensembl_gene_id, refseq_name).
#   file_mgi_annotations              — character, path to an RDS holding a GRanges of MGI TSS windows (metadata: feature_type, mgi_marker_id).
#   file_dhs_granges                  — character, path to an RDS holding a GRanges of DHS regions (metadata column: DHS).
#   file_refseq_select_gtf            — character, path to the RefSeq Select GTF, used to derive merged exonic regions.
#
# Returns:
#   list of data.frames, one per annotation source, each in the shape returned by summarise_overlap_DHS_with_X_ignore_direction() / _directional().
get_DHS_map_to_stuff = function(file_refseq_select_first_3_exons,
                                file_mgi_annotations,
                                file_dhs_granges,
                                file_refseq_select_gtf){

  res_list = list()
  ## refseq select
  res_list[[length(res_list)+1]] = get_GRange_of_mm39_genes_general(file_refseq = file_refseq_select_first_3_exons,
                                                                    output_mode = "TSS",
                                                                    postprocess_expand = 0) %>%
    summarise_overlap_DHS_with_X_directional(grange_x = .,  "refseq_select_gene_TSS", file_dhs_granges = file_dhs_granges)

  res_list[[length(res_list)+1]] = get_GRange_of_mm39_genes_general(file_refseq = file_refseq_select_first_3_exons,
                                                                    output_mode = "TTS",
                                                                    postprocess_expand = 0) %>%
    summarise_overlap_DHS_with_X_ignore_direction(grange_x = .,  "refseq_select_gene_TTS", file_dhs_granges = file_dhs_granges)

  res_list[[length(res_list)+1]] = get_GRange_of_mm39_genes_general(file_refseq = file_refseq_select_first_3_exons,
                                                                    output_mode = "full",
                                                                    postprocess_expand = 0) %>%
    summarise_overlap_DHS_with_X_ignore_direction(grange_x = .,  "refseq_select_gene_full", file_dhs_granges = file_dhs_granges)

  ## mgi noncoding

  grange_mgi = readRDS(file_mgi_annotations)
  granges_non_genic_elements = grange_mgi %>%
    subset(feature_type %in% c("lincRNA gene", "lncRNA gene", "miRNA gene",
                               "pseudogene", "snoRNA gene", "snRNA gene"))

  granges_non_genic_elements$gene_id = granges_non_genic_elements$mgi_marker_id
  granges_non_genic_elements$refseq_name = granges_non_genic_elements$mgi_marker_id
  res_list[[length(res_list)+1]] = summarise_overlap_DHS_with_X_directional(grange_x = granges_non_genic_elements,  "mgi_ncRNA", file_dhs_granges = file_dhs_granges)


  ## mgi enhancer
  granges_enhancers = grange_mgi %>%
    subset(feature_type %in% c("enhancer"))
  granges_enhancers$gene_id = granges_enhancers$mgi_marker_id
  granges_enhancers$refseq_name = granges_enhancers$mgi_marker_id
  res_list[[length(res_list)+1]] = summarise_overlap_DHS_with_X_ignore_direction(grange_x = granges_enhancers,  "mgi_enhancer", file_dhs_granges = file_dhs_granges)

  ## exonic regions
  ### refseq select exons
  grange_exon_refseq_select = get_non_overlap_exonic_region_based_on_refseq_select(file_refseq_select_gtf)
  mcols(grange_exon_refseq_select) = data.frame(gene_id = as.character(grange_exon_refseq_select),
                                                refseq_name = as.character(grange_exon_refseq_select))
  res_list[[length(res_list)+1]] = summarise_overlap_DHS_with_X_ignore_direction(grange_x = grange_exon_refseq_select,  "exon_refseq_select", file_dhs_granges = file_dhs_granges)

  return(res_list)
}


# ---- label_DHS_with_everything ----------------------------------------------
# Pivot the per-source overlap tables into one wide DHS table, twice: `annots_<type>_strand<p|m|b>` carries the comma-separated gene ids overlapping each DHS for each annotation type on each strand, and `filter_<type>_strand<p|m|b>` carries the mapping label for the same key, one of UNIQ, DHS-MULTI_ANNOT, ANNOT-MULTI_DHS, MULTI-MULTI, or NONE when nothing overlaps. The strand is folded to p/m/b (plus / minus / both), and the two pivots are joined on DHS.
#
# Parameters: the four paths taken by get_DHS_map_to_stuff(), which this function calls.
#
# Returns:
#   data.frame, one row per DHS: the annots_* and filter_* columns. A DHS with no overlap in a given source carries "NONE".
label_DHS_with_everything = function(file_refseq_select_first_3_exons,
                                     file_mgi_annotations,
                                     file_dhs_granges,
                                     file_refseq_select_gtf){

  list_overlap_DHS = get_DHS_map_to_stuff(file_refseq_select_first_3_exons = file_refseq_select_first_3_exons,
                                          file_mgi_annotations = file_mgi_annotations,
                                          file_dhs_granges = file_dhs_granges,
                                          file_refseq_select_gtf = file_refseq_select_gtf)

  df_annot_wide = list_overlap_DHS %>%
    bind_rows() %>%
    mutate(strand = case_when(
      strand == "+" ~ "p",
      strand == "-" ~ "m",
      strand == "*" ~ "b"
    ),
    annot_type_stranded = paste0(annot_type, "_strand", strand)) %>%
    group_by(DHS, annot_type_stranded) %>%
    summarise(n = n(),
              annot = paste0(gene_id, collapse = ","), .groups = "drop") %>%
    pivot_wider(id_cols = c("DHS"), names_from = "annot_type_stranded",
                values_from = annot,
                names_prefix = "annots_")

  df_annot_wide_filter = list_overlap_DHS %>%
    bind_rows() %>%
    mutate(strand = case_when(
      strand == "+" ~ "p",
      strand == "-" ~ "m",
      strand == "*" ~ "b"
    ),
    annot_type_stranded = paste0(annot_type, "_strand", strand)) %>%
    group_by(DHS, annot_type_stranded) %>%
    summarise(n = n(),
              annot_map_multi_DHS = any(annot_map_multi_DHS),
              DHS_map_multi_annot = any(DHS_map_multi_annot),
              annot = paste0(gene_id, collapse = ","), .groups = "drop") %>%
    mutate(label = case_when(
      annot_map_multi_DHS & DHS_map_multi_annot ~ "MULTI-MULTI",
      annot_map_multi_DHS  ~ "ANNOT-MULTI_DHS",
      DHS_map_multi_annot  ~ "DHS-MULTI_ANNOT",
      T ~ "UNIQ")) %>%
    pivot_wider(id_cols = "DHS", names_from = "annot_type_stranded",
                values_from = label,
                values_fill = "NONE",
                names_prefix = "filter_")

  stopifnot(all(df_annot_wide$DHS == df_annot_wide_filter$DHS))

  df_annot_merged = df_annot_wide %>%
    dplyr::inner_join(df_annot_wide_filter , by = "DHS")


  return(df_annot_merged)
}


# ---- str_remove_all_custom --------------------------------------------------
# Vectorised stringr::str_remove_all() that passes NA patterns through as the input string, so it can be applied to a column that mixes values and NA.
#
# Parameters:
#   string  — character vector to remove from.
#   pattern — character vector of patterns; NA leaves the input unchanged.
#
# Returns: character vector the same length as `string`.
str_remove_all_custom = Vectorize(function(string, pattern){

  ifelse(is.na(pattern), string,
         str_remove_all(string, pattern))

})


# ---- has_third_gene ---------------------------------------------------------
# Test whether an annotation string lists a gene other than the two genes named by the first- and second-strand TSS columns. Both names are stripped from the annotation string (all their occurrences), the remaining commas are removed, and the call returns TRUE when anything is left.
#
# Parameters:
#   first_strand_genename  — character vector, gene name(s) on the + strand.
#   second_strand_genename — character vector, gene name(s) on the - strand.
#   genenames_to_compare   — character vector, the annotation string to test, the same length as the two name vectors.
#
# Returns: logical vector the same length as the inputs.
has_third_gene = function(first_strand_genename, second_strand_genename,
                          genenames_to_compare) {
  stopifnot(length(first_strand_genename) ==
              length(second_strand_genename))
  stopifnot(length(first_strand_genename) == length(genenames_to_compare))

  # Remove first and second gene names from the list
  remaining = genenames_to_compare
  remaining = str_remove_all_custom(remaining, fixed(first_strand_genename))
  remaining = str_remove_all_custom(remaining, fixed(second_strand_genename))

  # Remove commas and check if anything remains
  remaining = str_remove_all_custom(remaining, ",")

  return(!is.na(remaining) & str_count(remaining) > 0)
}


# ---- curate_annotation_of_all_DHS -------------------------------------------
# Classify every DHS from the wide annotation table. It derives the exonic / intronic / intergenic position, whether a third gene's body overlaps, and then applies the ordered class rule:
#
#   annotation (full)          annotation_simp          condition
#   Known_1.PROMPT-Gene        1.Unique Gene            unique RefSeq-Select TSS, no third gene, no ncRNA
#   Known_2.Gene-Gene          2.Gene Filtered          any RefSeq-Select gene TSS
#   3.MGI-Enhancer             3.MGI-Enhancer           MGI enhancer AND non-exonic
#   4.OtherTSS                 4.OtherTSS               any other annotation
#   Unknown_6..8               5.Unknown                no known annotation
#
# Because annotation_simp is an ordered case_when, a lower-numbered class wins.
#
# Parameters:
#   file_dhs_annotation — character, path to an RDS holding the wide table written by label_DHS_with_everything(): one row per DHS, with the annots_* and filter_* columns.
#
# Returns:
#   data.frame, the wide annotation table plus annotation, annotation_simp, Known_TSS, Detail and Class columns. One row per DHS that overlaps at least one annotated element.
curate_annotation_of_all_DHS = function(file_dhs_annotation){
  df_annot_merged = readRDS(file_dhs_annotation)
  df_annot_merged_filtered = df_annot_merged %>%
    mutate(
      annotation_exon =
        filter_exon_refseq_select_strandb != "NONE",

      annotation_genebody =  filter_refseq_select_gene_full_strandb != "NONE" ,

      annotation_intronic = annotation_genebody & !annotation_exon,
      annotation_integenic = !annotation_genebody,

      annotation_third_gene = has_third_gene(annots_refseq_select_gene_TSS_strandp,
                                             annots_refseq_select_gene_TSS_strandm,
                                             annots_refseq_select_gene_full_strandb),


    )


  df_annot_merged_filtered_annotated = df_annot_merged_filtered %>%
    mutate(
      annotation_unique_gene_TSS =
        (filter_refseq_select_gene_TSS_strandp == "UNIQ" &
           filter_refseq_select_gene_TSS_strandm == "NONE" |
           filter_refseq_select_gene_TSS_strandm == "UNIQ" &
           filter_refseq_select_gene_TSS_strandp == "NONE") &
        # Not filtered on overlap with another gene's TTS.
        # More stringent: exclude overlap with another gene's body.
        ! annotation_third_gene &
        # dont overlap with ncRNA
        filter_mgi_ncRNA_strandp == "NONE" &
        filter_mgi_ncRNA_strandm == "NONE",  # DHS not map to ncRNA

      annotation_stranded_gene_TSS =
        (filter_refseq_select_gene_TSS_strandp == "UNIQ" |
           filter_refseq_select_gene_TSS_strandm == "UNIQ") &

        # dont overlap with other TTS
        # Not filtered on overlap with another gene's TTS.
        # dont overlap with gene body of another gene
        !annotation_third_gene &

        filter_mgi_ncRNA_strandp == "NONE" &
        filter_mgi_ncRNA_strandm == "NONE",

      annotation_any_refseq_select_gene_TSS =
        (filter_refseq_select_gene_TSS_strandp != "NONE" |
           filter_refseq_select_gene_TSS_strandm != "NONE") ,


      annotation_MGIenhancer =
        (annotation_intronic | annotation_integenic) &
        (filter_mgi_enhancer_strandb != "NONE"),

      annotation_any_annotation =
        filter_refseq_select_gene_TSS_strandm != "NONE" |
        filter_refseq_select_gene_TSS_strandp != "NONE" |
        filter_mgi_ncRNA_strandm != "NONE" |
        filter_mgi_ncRNA_strandp != "NONE" |
        filter_mgi_enhancer_strandb != "NONE",

      annotation_unknown_exonic =
        !annotation_any_annotation &
        annotation_exon,

      annotation_unknown_intronic = (
        !annotation_any_annotation &
          annotation_intronic),

      annotation_unknown_integenic = (
        !annotation_any_annotation &
          annotation_integenic)
    ) %>%
    mutate(annotation =
             case_when(
               annotation_unique_gene_TSS ~ "Known_1.PROMPT-Gene",
               annotation_stranded_gene_TSS ~ "Known_2.Gene-Gene",

               annotation_any_refseq_select_gene_TSS ~ "Known_3.GeneFiltered",
               annotation_MGIenhancer ~ "Known_4.MGI-Enhancer",
               annotation_any_annotation ~ "Known_5.OtherTSS",
               annotation_unknown_exonic ~ "Unknown_6.Exonic",
               annotation_unknown_intronic ~ "Unknown_7.Intronic",
               annotation_unknown_integenic ~ "Unknown_8.Integenic",
               T ~ "Unknown_Unknown"
             ),
           annotation_simp = case_when(
             # The enhancer class below is the plain MGI-enhancer rule and does not test annotation_third_gene; that flag gates only the two gene classes above.
             annotation_unique_gene_TSS ~ "1.Unique Gene",
             annotation_any_refseq_select_gene_TSS | annotation_stranded_gene_TSS ~ "2.Gene Filtered",
             annotation_MGIenhancer ~ "3.MGI-Enhancer",
             annotation_any_annotation ~ "4.OtherTSS",
             annotation_unknown_exonic | annotation_unknown_intronic | annotation_unknown_integenic ~ "5.Unknown",
             T ~ "Unknown_Unknown"
           )
    ) %>%
    mutate(Known_TSS = str_extract(annotation, "(\\w+)_([\\w-\\._]+)", 1),
           Detail = str_extract(annotation, "(\\w+)_([\\w-\\._]+)", 2),
           Class = case_when(
             str_detect(annotation, "Gene") ~ "Genic",
             str_detect(annotation, "Enhancer") ~ "Enhancer",
             str_detect(annotation, "Unknown") ~ "Unknown",
             T ~ "Others"
           ),
           Class = factor(Class, levels = c("Genic", "Enhancer", "Others", "Unknown"))) %>%
    mutate(Class = factor(Class))


  return(df_annot_merged_filtered_annotated)
}


# ---- get_DHS_map_to_enhancer ------------------------------------------------
# Select the enhancer class from the classified annotation table.
#
# Parameters:
#   file_dhs_classification — character, path to an RDS holding the table returned by curate_annotation_of_all_DHS(), with its annotation_simp column.
#
# Returns:
#   data.frame, the rows of the classified table whose annotation_simp is "3.MGI-Enhancer".
get_DHS_map_to_enhancer = function(file_dhs_classification){
  readRDS(file_dhs_classification) %>%
    filter(annotation_simp == "3.MGI-Enhancer")
}

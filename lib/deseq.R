## DESeq2 with size factors computed from chrM read counts.
## Adapted from the DHS annotation code (03-Annotation).
## Requires CHRM_DEPTH (config/paths.R), tmp_replicates (lib/tldr_reads.R) and
## define_scale_factor() (04-Normalization/scale_factors.R).

library(tidyverse)
library(DESeq2)

#' run deseq on count matrix with scale factor calculated from mitochondrial reads
#' @param df_count_matrix dataframe, where the columns are sample ids, e.g. mES_WT_rep1
#' @param df_meta dataframe, where rownames are sample ids (same order as the count matrix columns)
#' @return list: df_lfc_shrunk (apeglm), dds, df_lfc_raw; one row per feature and genotype vs WT
run_deseq <- function(df_count_matrix, df_meta,
                      candidate_geno = c("ARS2", "RBM7", "ZC3H18", "ZCCHC8"),
                      df_scale_factor = define_scale_factor(CHRM_DEPTH, sample_levels = tmp_replicates)) {
    stopifnot(all(colnames(df_count_matrix) == rownames(df_meta)))
    stopifnot(all(str_detect(colnames(df_count_matrix), "^mES_([A-Z0-9]+)_(rep\\d)$")))
    stopifnot("geno" %in% colnames(df_meta))
    stopifnot(is.factor(df_meta$geno))
    stopifnot(levels(df_meta$geno)[1] == "WT")
    stopifnot(all(candidate_geno %in% levels(df_meta$geno)))

    if (!is.null(df_scale_factor)) {
        df_scale_factor_ordered <- df_meta %>%
            rownames_to_column("sample_id") %>%
            dplyr::left_join(df_scale_factor, by = c("sample_id" = "sample")) %>%
            dplyr::select(c("sample_id", scale_factor))
        stopifnot(all(colnames(df_count_matrix) == df_scale_factor_ordered$sample_id))
    }

    dds <- DESeqDataSetFromMatrix(countData = df_count_matrix,
                                  colData = df_meta,
                                  design = ~ geno)
    if (!is.null(df_scale_factor)) {
        sizeFactors(dds) <- df_scale_factor_ordered$scale_factor
    }
    dds <- estimateDispersions(dds)
    dds <- nbinomWaldTest(dds)

    tmp_list_raw <- list()
    for (i in candidate_geno) {
        tmp_list_raw[[i]] <- DESeq2::results(dds,
                                             contrast = c("geno", i, "WT"),
                                             pAdjustMethod = "fdr",
                                             lfcThreshold = 0) %>%
            as.data.frame() %>%
            rownames_to_column("DHS") %>%
            mutate(metric = paste0("lfc_", i))
    }
    df_lfc_raw <- bind_rows(tmp_list_raw)

    tmp_list_shrunk <- list()
    for (i in candidate_geno) {
        tmp_list_shrunk[[i]] <- DESeq2::lfcShrink(dds,
                                                  coef = paste0("geno_", i, "_vs_WT"),
                                                  type = "apeglm") %>%
            as.data.frame() %>%
            rownames_to_column("DHS") %>%
            mutate(metric = paste0("lfc_", i))
    }
    df_lfc_shrunk <- bind_rows(tmp_list_shrunk)

    return(list(df_lfc_shrunk = df_lfc_shrunk,
                dds = dds,
                df_lfc_raw = df_lfc_raw))
}

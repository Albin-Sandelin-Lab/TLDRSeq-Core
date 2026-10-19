## Read the TLDR-Seq reads assigned to each DHS.
## Adapted from the DHS annotation code (03-Annotation).
## Requires DIR_READS_PER_DHS (config/paths.R).
##
## Each read is one line of <DIR_READS_PER_DHS>/<DHS>/<replicate>/reads_score_by_pA_length.bed.
## Columns used downstream: V1 chrom, V2 start, V3 end, V5 pA tail length, V6 strand,
## V10 number of exons, V11 exon lengths, V13 DHS, V14 replicate.

library(tidyverse)
library(pbapply)

## TLDRSeq_sup samples
tmp_replicates <- c("mES_WT_rep1", "mES_WT_rep2", "mES_WT_rep3", "mES_WT_rep4", "mES_WT_rep5", "mES_WT_rep6",
                    "mES_ARS2_rep1", "mES_ARS2_rep2", "mES_ARS2_rep3",
                    "mES_ZC3H18_rep1", "mES_ZC3H18_rep2", "mES_ZC3H18_rep3",
                    "mES_ZC3H4_rep1", "mES_ZC3H4_rep2", "mES_ZC3H4_rep3",
                    "mES_RBM7_rep1", "mES_RBM7_rep2", "mES_RBM7_rep3",
                    "mES_ZCCHC8_rep1", "mES_ZCCHC8_rep2", "mES_ZCCHC8_rep3")

## read.table that returns an empty data.frame for empty files
safe_read_table <- function(file, ...) {
    tryCatch(
    {
        data <- read.table(file, ...)
        return(data)
    }, error = function(e) {
        if (grepl("no lines available in input", e$message)) {
            return(data.frame())
        } else {
            stop(e)
        }
    })
}

read_TLDR_reads_of_one_DHS_in_one_sample <- function(query_DHS, query_replicate, dir_files = DIR_READS_PER_DHS) {
    stopifnot(length(query_DHS) == 1)
    stopifnot(length(query_replicate) == 1)
    tmp_dir_dhs <- file.path(dir_files, query_DHS, query_replicate)
    tmp_file_reads <- file.path(tmp_dir_dhs, "reads_score_by_pA_length.bed")
    if (!file.exists(tmp_file_reads)) {
        ## "empty" marker = DHS not expressed in this sample (expected)
        if (file.exists(file.path(tmp_dir_dhs, "empty"))) {
            return(data.frame())
        } else {
            stop(sprintf("DHS %s is not quantified in sample %s. Neither BED file nor 'empty' marker found in: %s",
                         query_DHS, query_replicate, tmp_dir_dhs))
        }
    } else {
        return(safe_read_table(tmp_file_reads,
                               colClasses = c("character", "integer", "integer", "character",
                                              "integer", "character", "integer", "integer",
                                              "character", "integer", "character",
                                              "character", "character", "character")))
    }
}

read_TLDR_reads_by_list_of_DHS <- function(query_list_DHS, query_list_replicates, dir_files = DIR_READS_PER_DHS) {
    combinations <- expand.grid(DHS = unique(query_list_DHS), replicate = unique(query_list_replicates))
    arg_list <- split(combinations, seq(nrow(combinations)))
    tmp_list <- pblapply(arg_list, function(args) {
        read_TLDR_reads_of_one_DHS_in_one_sample(args$DHS, args$replicate, dir_files = dir_files)
    })
    return(bind_rows(tmp_list))
}

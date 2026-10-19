## Label the 3' end of the TLDRSeq_sup reads of each gene-associated DHS, with
## label_reads_by_3prime_categories() (label_reads_by_3prime.R), and count the reads of each
## label combination per DHS and replicate.
##
## Usage (from 06-Genic/read-classes/): Rscript label_reads_TLDRSeq_sup.R [n_cores]
##   n_cores : DHSs processed in parallel (default 8)
##
## Input : DHS_GENE_MAP_EXONS, DHS_TSS_CURATED, DHS_TTS_CURATED, DIR_READS_PER_DHS (config/paths.R)
## Output: TTS_PROP_CATEGORIES_RDS (config/paths.R)
##         one row per DHS x replicate x label combination: labels (label_tts, label_cdsts,
##         label_premature_termination, unspliced, label_premature, label_newly_degraded,
##         label_tts_falling_region, polyA_tail), reference positions of the gene (typical_*),
##         read summaries (mean / median length, exon count, pA length, 5' and 3' ends,
##         distance to the reference TSS and TTS), n (reads) and freq (n / reads of the DHS
##         in the replicate)
##
## Windows around the reference TTS and CDS end: +/- 50 bp. pA+: poly(A) tail > 15 nt.

library(parallel)
library(GenomicRanges)
library(pbapply)
library(tidyverse)

source("../../config/paths.R")
source("../../lib/tldr_reads.R")                 # tmp_replicates: the 21 TLDRSeq_sup replicates
source("../../02-DHS/raw_reads_access.R")
source("exon_structure.R")
source("label_reads_by_3prime.R")

args <- commandArgs(trailingOnly = TRUE)
n_cores <- ifelse(length(args) >= 1, as.integer(args[1]), 8)

dir.create(dirname(TTS_PROP_CATEGORIES_RDS), showWarnings = FALSE, recursive = TRUE)


#### Gene features of each DHS

df_DHS_gene_unique_map_with_exon_info <- readRDS(DHS_GENE_MAP_EXONS)
grange_tts <- readRDS(DHS_TTS_CURATED)
grange_tss <- readRDS(DHS_TSS_CURATED)

## CDS end and start of the gene, in the gene orientation
grange_cdsts <- df_DHS_gene_unique_map_with_exon_info %>%
    transmute(seqnames = chrom,
              start = ifelse(strand == "+", cdsEnd, cdsStart),
              end = ifelse(strand == "+", cdsEnd, cdsStart),
              strand = strand,
              gene_id = entrezgene_id,
              gene_name,
              refseq_name,
              DHS) %>%
    GRanges()

grange_cdsss <- df_DHS_gene_unique_map_with_exon_info %>%
    transmute(seqnames = chrom,
              start = ifelse(strand == "+", cdsStart, cdsEnd),
              end = ifelse(strand == "+", cdsStart, cdsEnd),
              strand = strand,
              gene_id = entrezgene_id,
              gene_name,
              refseq_name,
              DHS) %>%
    GRanges()

stopifnot(!any(is.na(grange_tts$gene_name)))
stopifnot(!any(is.na(grange_cdsts$gene_name)))
stopifnot(all(grange_tts$gene_name %in% grange_cdsts$gene_name))
stopifnot(all(grange_cdsts$gene_name %in% grange_tts$gene_name))


#### Label the reads of each DHS x replicate, count per label combination

label_reads_of_one_DHS <- function(i) {
    tmp_gene_name <- df_DHS_gene_unique_map_with_exon_info %>%
        filter(DHS == i) %>%
        .$gene_name

    if (!(length(tmp_gene_name) == 1)) {
        return(NULL)                             # DHS mapped to several genes
    }
    if (!(tmp_gene_name %in% grange_tts$gene_name)) {
        return(NULL)                             # no reference TTS
    }
    if (!(tmp_gene_name %in% grange_cdsts$gene_name)) {
        return(NULL)                             # no CDS end
    }

    tmp_typical_tss <- grange_tss[grange_tss$DHS == i] %>% start()
    tmp_typical_tts <- grange_tts[grange_tts$DHS == i] %>% start()
    tmp_typical_cdsss <- grange_cdsss[grange_cdsss$DHS == i] %>% start()
    tmp_typical_cdsts <- grange_cdsts[grange_cdsts$DHS == i] %>% start()

    stopifnot(is.numeric(tmp_typical_tss))
    stopifnot(is.numeric(tmp_typical_tts))
    stopifnot(is.numeric(tmp_typical_cdsss))
    stopifnot(is.numeric(tmp_typical_cdsts))

    tmp_list <- list()
    for (j in tmp_replicates) {
        tmp_processed_data <- label_reads_by_3prime_categories(query_DHS = i,
                                                               query_sample = j,
                                                               grange_tts = grange_tts,
                                                               grange_cds = grange_cdsts,
                                                               df_DHS_gene_map_with_exon_info = df_DHS_gene_unique_map_with_exon_info,
                                                               dir_dhs_depth = DIR_READS_PER_DHS,
                                                               postion_window_size = 100)
        if (nrow(tmp_processed_data) == 0) {
            next
        }
        tmp_list[[length(tmp_list) + 1]] <- tmp_processed_data %>%
            mutate(sample = j) %>%
            group_by(gene_name,
                     DHS,
                     sample,
                     strand,
                     label_tts,
                     label_cdsts,
                     label_premature_termination,
                     unspliced = exon_count == 1,
                     label_premature = label_unspliced_premature_termination == "premature_unspliced",
                     label_newly_degraded,
                     label_tts_falling_region,
                     polyA_tail = polyA_length > 15) %>%
            summarise(typical_exon_count = median(typical_exon_count),
                      typical_tss = tmp_typical_tss,
                      typical_tts = tmp_typical_tts,
                      typical_cdsss = tmp_typical_cdsss,
                      typical_cdsts = tmp_typical_cdsts,
                      mean_exonLength = mean(exonLength),
                      median_exonLength = median(exonLength),
                      mean_length = mean(length),
                      median_length = median(length),
                      mean_exon_count = mean(exon_count),
                      median_exon_count = median(exon_count),
                      median_polyA_tail_length = median(polyA_length),
                      mean_polyA_tail_length = mean(polyA_length),
                      median_TTS = median(TTS),
                      mean_TTS = mean(TTS),
                      median_TSS = median(TSS),
                      mean_TSS = mean(TSS),
                      median_distance_to_TTS = median(TTS - typical_tts),
                      mean_distance_to_TTS = mean(TTS - typical_tts),
                      median_distance_to_TSS = median(TSS - typical_tss),
                      mean_distance_to_TSS = mean(TSS - typical_tss),
                      n = n(),
                      .groups = "drop") %>%
            mutate(freq = n / sum(n))
    }
    bind_rows(tmp_list)
}

list_per_DHS <- unique(df_DHS_gene_unique_map_with_exon_info$DHS) %>%
    mclapply(label_reads_of_one_DHS, mc.cores = n_cores)
stopifnot(!any(sapply(list_per_DHS, inherits, "try-error")))

df_TTS_prop_categories <- bind_rows(list_per_DHS)

saveRDS(df_TTS_prop_categories, TTS_PROP_CATEGORIES_RDS)

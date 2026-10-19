## Library depth and chrM scale factor of each replicate, from bambu gene counts.
##
## Usage (from 06-Genic/read-classes/): Rscript 02_library_size_factors.r <dataset> [n_cores]
##   dataset : TLDRSeq_nuclear_sup (library depth, used for CPM) or
##             TLDRSeq_Z1_MTR4_sup (chrM scale factors, used by DESeq2)
##   n_cores : bambu cores (default: number of BAM files)
##
## Input : $RESULTS_DIR/<dataset>/concatenate_per_replicate/<sample>.bam (01-TLDRSeq-alignment)
## Output (in $GENIC_DIR/<dataset>/library_size_factors/):
##   gene_exp.tsv               bambu gene counts and CPM per sample
##   library_size_factors.tsv   sample, depth (sum of gene counts), chrM_counts (gene counts on
##                              chrMT), scale_factor_chrM (chrM_counts / mean chrM_counts)

library(tidyverse)
library(rtracklayer)
library(bambu)

source("../../config/paths.R")

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]

in_dir  <- file.path(RESULTS_DIR, dataset, "concatenate_per_replicate")
out_dir <- file.path(GENIC_DIR, dataset, "library_size_factors")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

bam.name <- file.path(in_dir, list.files(in_dir, pattern = "\\.bam$"))
n_cores <- ifelse(length(args) >= 2, as.integer(args[2]), length(bam.name))


#### Annotation: the genome names the mitochondrial chromosome chrMT

gtf <- rtracklayer::import(GENCODE_GTF) %>%
    renameSeqlevels(c("chrM" = "chrMT"))
gtf.adjusted.name <- file.path(out_dir, "fixed_gencode_annotation.gtf")
export(gtf, gtf.adjusted.name, format = "gtf")

gtf.gene.names <- gtf %>% as.data.frame %>%
    dplyr::select(gene.id = gene_id, gene.name.gtf = gene_name, seqnames) %>%
    distinct %>%
    mutate(gene.id = gene.id %>% str_replace("\\..*", ""))


#### bambu gene counts

annotations <- prepareAnnotations(gtf.adjusted.name)

se.quant <- bambu(reads = bam.name, annotations = annotations,
                  genome = GENOME_FA, discovery = FALSE, ncore = n_cores)

gse.quant <- se.quant %>%
    transcriptToGeneExpression()

gse.df <- bind_rows(CPM = reshape2::melt(assays(gse.quant)$CPM),
                    counts = reshape2::melt(assays(gse.quant)$counts),
                    .id = "metric")
write.table(gse.df, file = file.path(out_dir, "gene_exp.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


#### Library depth and chrM scale factors

depth <- gse.df %>%
    dplyr::filter(metric == "counts") %>%
    group_by(sample = Var2) %>%
    dplyr::summarize(depth = sum(value))

mito.counts <- gse.df %>%
    mutate(gene.id = Var1 %>% str_replace("\\..*", "")) %>%
    dplyr::filter(metric == "counts") %>%
    left_join(gtf.gene.names, by = "gene.id") %>%
    dplyr::filter(seqnames == "chrMT") %>%
    group_by(sample = Var2) %>%
    dplyr::summarize(chrM_counts = sum(value)) %>%
    mutate(scale_factor_chrM = chrM_counts / mean(chrM_counts))

size.factors <- depth %>%
    left_join(mito.counts, by = "sample")
write.table(size.factors, file = file.path(out_dir, "library_size_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

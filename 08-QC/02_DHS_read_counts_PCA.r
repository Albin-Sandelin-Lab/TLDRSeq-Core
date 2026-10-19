## PCA of the TLDR-Seq read counts per DHS and strand (TLDRSeq_sup: WT, RBM7, ZCCHC8).
##
## Usage (from 08-QC/): Rscript 02_DHS_read_counts_PCA.r [n_cores]
##   n_cores : DHSs read in parallel (default 15)
##
## All DHSs of the annotation; reads counted per replicate, DHS and strand. The 2 DHSs with
## the most reads are removed. PCA on log(1 + count), centered and scaled.
##
## Input : DHS_ANNOTATION_RDS, DIR_READS_PER_DHS (config/paths.R)
## Output (in $QC_DIR/DHS_read_counts_PCA/):
##   DHS_reads_counts_nopA.tsv        rep, DHS, strand, N (all replicates of TLDRSeq_sup)
##   DHS_read_counts_PC1_2_plot.pdf   PC1 vs PC2 (WT, RBM7, ZCCHC8)

library(parallel)
library(data.table)
library(tidyverse)

source("../config/paths.R")
source("../lib/tldr_reads.R")

args <- commandArgs(trailingOnly = TRUE)
n_cores <- ifelse(length(args) >= 1, as.integer(args[1]), 15)

out_dir <- file.path(QC_DIR, "DHS_read_counts_PCA")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#### Reads per replicate, DHS and strand

dhs.list <- readRDS(DHS_ANNOTATION_RDS)$DHS %>% unique

count.dhs.reads <- function(dhs) {
    DHS_reads <- read_TLDR_reads_by_list_of_DHS(query_list_DHS = dhs, query_list_replicates = tmp_replicates)
    if (nrow(DHS_reads) == 0) {
        return(NULL)
    }
    as.data.table(DHS_reads)[, .N, by = .(rep = V14, DHS = V13, strand = V6)]
}

dhs.counts.list <- mclapply(dhs.list, count.dhs.reads, mc.cores = n_cores)
dhs.reads.counts.nopA <- rbindlist(dhs.counts.list[!sapply(dhs.counts.list, is.null)])
rm(dhs.counts.list)

## remove the 2 DHSs with the most reads
remove.dhs <- dhs.reads.counts.nopA[, .(n = sum(N)), by = DHS][order(-n)]$DHS %>% head(2)
dhs.reads.counts.nopA <- dhs.reads.counts.nopA[!DHS %in% remove.dhs, ]

fwrite(dhs.reads.counts.nopA, file.path(out_dir, "DHS_reads_counts_nopA.tsv"), sep = "\t")


#### PCA

dhs.prep.FC <- dhs.reads.counts.nopA %>%
    mutate(DHS.strand = paste(DHS, strand, sep = ".")) %>%
    dplyr::select(DHS.strand, n = N, rep) %>%
    dplyr::filter(!grepl("ZC3H4|ZC3H18|ARS2", rep)) %>%
    reshape2::dcast(DHS.strand ~ rep, value.var = "n", fill = 0) %>%
    column_to_rownames("DHS.strand")

dhs.prep.FC.mat <- as.matrix(dhs.prep.FC) %>%
    {log(1 + .)}
pca.obj <- prcomp(dhs.prep.FC.mat %>% t, scale. = TRUE, center = TRUE)

axis_lab <- paste(as.character(colnames(summary(pca.obj)$importance)), " (", as.character(round(summary(pca.obj)$importance[2, ] * 100)), "%)", sep = "")
tab <- as.data.frame(pca.obj$x) %>%
    rownames_to_column("sampleID") %>%
    mutate(geno = sampleID %>% str_replace("mES_", "") %>% str_replace("_rep.", "")) %>%
    mutate(geno = factor(geno, levels = c("WT", "RBM7", "ZCCHC8")))

p1 <- ggplot(tab, aes(x = PC1, y = PC2, colour = geno)) +
    geom_point(size = 3) +
    theme_bw() +
    scale_colour_manual(values = c('#828282', "#3D405B", "#F2CC8F")) +
    xlab(axis_lab[1]) +
    ylab(axis_lab[2]) +
    ggtitle("PCA of DHS read counts")
ggsave(p1, filename = file.path(out_dir, "DHS_read_counts_PC1_2_plot.pdf"), width = 6, height = 4)

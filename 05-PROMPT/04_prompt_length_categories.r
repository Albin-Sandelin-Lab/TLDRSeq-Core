## PROMPT sensitivity to RBM7 / ZCCHC8 depletion by read length category.
##
## Usage (from 05-PROMPT/): Rscript 04_prompt_length_categories.r
##
## DESeq2 (WT, RBM7, ZCCHC8; chrM size factors) on PROMPT, PMT, exon1 and gene read counts
## per DHS and genomic read length category; plots for PROMPT reads:
##   - fraction of upregulated PROMPTs (log2FC > 1, padj < 0.05) per length category, pA+ / pA-
##     (DHSs with > 2 normalized pA+ or pA- reads in the category)
##   - log2FC per length category (violin plots)
##
## Input : $PROMPT_DIR/prompt_vs_pmt_length/ (03_prompt_vs_pmt_length.r)
## Output: $PROMPT_DIR/prompt_length_categories/
##   DHS_up_fraction_barplot_length_separated_pA_nopA_prompt.pdf
##   DHS_log2FC_violin_plot_length_separated_prompt.pdf
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(parallel)
library(pheatmap)
library(data.table)
library(org.Mm.eg.db)
library(tidyverse)
library(RColorBrewer)
library(classInt)
library(ggVennDiagram)

source("../config/paths.R")
source("../lib/tldr_reads.R")
source("../lib/deseq.R")
source("../04-Normalization/scale_factors.R")

out_dir <- file.path(PROMPT_DIR, "prompt_length_categories")
dir.create(out_dir, showWarnings=FALSE, recursive =TRUE)

out_dir.bk <- file.path(PROMPT_DIR, "prompt_vs_pmt_length")
dhs.reads.number.summarized.pA <- fread( file.path(out_dir.bk, "DHS_reads_number_summarized_pA.tsv"))
dhs.reads.statistics.pA <- fread( file.path(out_dir.bk, "DHS_reads_statistics_pA.tsv"))
dhs.reads.number.summarized.merged <- fread( file.path(out_dir.bk, "DHS_reads_number_summarized_merged.tsv"))
dhs.reads.statistics.merged <- fread( file.path(out_dir.bk, "DHS_reads_statistics_merged.tsv"))
dhs.reads.numbers.nopA <- fread(file.path(out_dir.bk, "DHS_reads_counts_nopA.tsv"))

#### sensitivity related
#### prepare for DESeq2

dhs.prep.FC <- dhs.reads.numbers.nopA %>%
    dplyr::select(dhs.type, group_cat, n = N, rep) %>%
    mutate(dhs.type.group_cat = paste0(dhs.type,".", group_cat)) %>% 
    dplyr::filter(!grepl("ZC3H4|ZC3H18|ARS2", rep)) %>% 
    reshape2::dcast(dhs.type.group_cat~ rep, value.var = "n", fill = 0) %>% 
    column_to_rownames("dhs.type.group_cat") ## %>%
    ## setNames(paste0("scale_",colnames(.)))

df.meta <- colnames(dhs.prep.FC) %>%
    data.frame(sample = .) %>%
    mutate(geno = sample %>% str_replace("mES_","") %>% str_replace("_rep.",""))  %>%
    column_to_rownames("sample") %>%
    mutate(geno = factor(geno, c("WT",  "RBM7",  "ZCCHC8")))

deseq.results <- run_deseq(dhs.prep.FC, df.meta, candidate_geno =c("RBM7",  "ZCCHC8"))

deseq.results.df <- deseq.results$df_lfc_shrunk  %>%
    dplyr::select(DHS, log2FoldChange , metric, padj) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    mutate(geno = factor(geno, levels = c("RBM7",   "ZCCHC8")))

deseq.results.df_up <- deseq.results.df %>%
    drop_na %>% 
    dplyr::filter(log2FoldChange > 1, padj < 0.05 ) %>%
    dplyr::select(DHS, geno) %>%
    mutate(status = "upredulated") %>%
    separate(DHS, into = c("DHS", "type", "group_cat"), sep = "\\.")

group_cat.levels <- c("0-100", "100-150", "150-200", "200-250", "250-300", "300-350", "350-400",
                      "400-500", "500-600", "600-700",  "700-800" , "800-1000", "1000-1500",  ">1500")

deseq.results.df.joined <- deseq.results.df %>%
    separate(DHS , into = c("DHS", "type", "group_cat"), sep = "\\.") %>%
    left_join(dhs.reads.number.summarized.merged %>%
              separate(dhs.type, into = c("DHS","type"), sep="\\."),
              by = c("DHS", "geno","type", "group_cat")) %>%
    dplyr::filter(!is.na(scaled_total_per_lib))  %>%
    left_join(dhs.reads.number.summarized.pA %>%
              separate(dhs.type, into = c("DHS","type"), sep="\\.") %>%
              dplyr::select(-strand, -total_per_lib, -total, -scaled_total),
              by = c("DHS", "geno", "group_cat", "type"),
              multiple ="all") %>%
    dplyr::rename( tpm.merged = scaled_total_per_lib.x,
                  tpm.pA = scaled_total_per_lib.y) %>%
    dplyr::filter(tpm.pA > 2) %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>% 
    mutate(status = ifelse(padj < 0.05 & log2FoldChange > 1, "upregulated", "unchanged"))  

### prepare to plot
deseq.results.df.joined.summarized <- deseq.results.df.joined %>%
    group_by(geno, group_cat, status, pA_status, type) %>%
    dplyr::summarize(  n = n() ) %>%
    group_by(geno, group_cat, pA_status, type) %>%
    mutate( n.total = sum(n), fraction = n/n.total) %>%
    mutate(group_cat =
               factor(group_cat,
                      levels = group_cat.levels)) %>%
    mutate(geno = factor(geno, levels = c( "RBM7",  "ZCCHC8")))  

## for (elt in c("pmt", "prompt", "gene", "exon1")) {
for (elt in c("prompt")) {
    deseq.results.df.joined.summarized.elt <- deseq.results.df.joined.summarized %>%
        dplyr::filter(type == elt)
    g <- ggplot(deseq.results.df.joined.summarized.elt %>%
              dplyr::filter(status =="upregulated"), aes(x = geno, y = fraction, fill = geno)) +
        geom_bar(stat ="identity" , color = "black", size = 0.2) + 
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x  = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 6),
              strip.placement = "outside") +
        facet_grid(vars( pA_status), vars(group_cat), switch = "x") +
        scale_fill_manual(values = c("#3D405B","#F2CC8F")) +
        xlab("DHS length category") +
        ylab("fraction of upregulated DHS") +
        ggtitle(paste0("DHS upregulated count upon exosome related factor KD - ", elt, " reads"))
    ggsave(g, filename = file.path(out_dir, paste0("DHS_up_fraction_barplot_length_separated_pA_nopA_", elt, ".pdf")), width = 10, height =4)
}

### make violin plots
deseq.results.df.joined.plot.violin <-  deseq.results.df.joined %>%
    dplyr::select(-status, -tpm.pA, -tpm.merged, -pA_status) %>%
    distinct() %>%
    mutate(group_cat =
               factor(group_cat,
                      levels = group_cat.levels)) %>%
    mutate(geno = factor(geno, levels = c( "RBM7", "ZCCHC8")))

## for (elt in c("pmt", "prompt", "gene", "exon1")) {
for (elt in c("prompt")) {
    deseq.results.df.joined.plot.violin.elt <- deseq.results.df.joined.plot.violin %>%
        dplyr::filter(type == elt)
    g <- ggplot(deseq.results.df.joined.plot.violin.elt, aes( x = group_cat, y= log2FoldChange, fill = geno)) +
        geom_violin( color = "black", size = 0.2) +
        theme_bw() +
        theme(# axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x  = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 6),
            strip.placement = "outside") +
        stat_summary(fun =mean, geom="point", shape=1, size=2, color="black") +
        stat_summary(fun = median, geom="point", shape=2, size=2, color="black") +
        facet_grid(row = vars(geno),  switch = "x") +
        scale_fill_manual(values = c("#3D405B","#F2CC8F")) +
        xlab("DHS length category") +
        ylab("log2 Fold Change vs WT") +
        ggtitle(paste0("DHS log2FC upon exosome related factor KD - ", elt, " reads"))
    ggsave(g, filename = file.path(out_dir, paste0("DHS_log2FC_violin_plot_length_separated_", elt, ".pdf")), width = 10, height =4)
}


## PROMPT sensitivity to RBM7 / ZCCHC8 depletion.
##
## Usage (from 05-PROMPT/): Rscript 02_prompt_sensitivity.r
##
## Input : $PROMPT_DIR/prompt_reads/ (01_prompt_reads_and_stats.r)
## Output: $PROMPT_DIR/sensitivity/
##   deseq_results_DHS.rds    DESeq2 PROMPT log2FC (RBM7, ZCCHC8 vs WT, chrM size factors)
##   DHS_list_expressed_reads_over_5.txt
##   PDF figures: log2FC density, Venn of upregulated PROMPTs, RBM7 vs ZCCHC8 log2FC,
##   read length / exon number / pA length of upregulated PROMPTs, pA fraction,
##   pA fraction vs log2FC violin plots, sensitivity and exonic length heatmaps
##
## Upregulated PROMPT: log2FC > 1 and padj < 0.05. "Expressed": > 5 normalized reads.
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
library(ggExtra)


source("../config/paths.R")
source("../lib/tldr_reads.R")
source("../lib/deseq.R")
source("../04-Normalization/scale_factors.R")

out_dir <- file.path(PROMPT_DIR, "sensitivity")
dir.create(out_dir, showWarnings=FALSE, recursive =TRUE)

out_dir.bk <- file.path(PROMPT_DIR, "prompt_reads")
dhs.reads.number.summarized.pA <- fread( file.path(out_dir.bk, "DHS_reads_number_summarized_pA.tsv"))
dhs.reads.statistics.pA <- fread( file.path(out_dir.bk, "DHS_reads_statistics_pA.tsv"))
dhs.reads.number.summarized.merged <- fread( file.path(out_dir.bk, "DHS_reads_number_summarized_merged.tsv"))
dhs.reads.statistics.merged <- fread( file.path(out_dir.bk, "DHS_reads_statistics_merged.tsv"))
dhs.reads.numbers.nopA <- fread(file.path(out_dir.bk, "DHS_reads_counts_nopA.tsv"))

#### sensitivity related
#### prepare for DESeq2
dhs.prep.FC <- dhs.reads.numbers.nopA %>%
    mutate(DHS.strand = paste(DHS, strand, sep = ".")) %>%
    dplyr::select(DHS.strand, n = N, rep) %>%
    dplyr::filter(!grepl("ZC3H4|ZC3H18|ARS2", rep)) %>% 
    reshape2::dcast(DHS.strand ~ rep, value.var = "n", fill = 0) %>%
    column_to_rownames("DHS.strand") ## %>%
    ## setNames(paste0("scale_",colnames(.)))

df.meta <- colnames(dhs.prep.FC) %>%
    data.frame(sample = .) %>%
    mutate(geno = sample %>% str_replace("mES_","") %>% str_replace("_rep.",""))  %>%
    column_to_rownames("sample") %>%
    mutate(geno = factor(geno, c("WT",  "RBM7",  "ZCCHC8")))

deseq.results <- run_deseq(dhs.prep.FC, df.meta, candidate_geno = c(  "RBM7",  "ZCCHC8"))

saveRDS(deseq.results, file.path(out_dir, "deseq_results_DHS.rds"))

list.dhs.tested <- deseq.results$df_lfc_shrunk$DHS %>% unique

## plot log2FC 
deseq.results.df <- deseq.results$df_lfc_shrunk  %>%
    dplyr::select(DHS, log2FoldChange , metric, padj) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    mutate(geno = factor(geno))

### make same plot for highly >5 reads prompts
deseq.results.df_for.density <- deseq.results.df %>%
    dplyr::filter(geno != "ZC3H4") %>%
    mutate(strand = str_replace(DHS, ".*\\.", ""))  %>%
    mutate(DHS = str_replace(DHS, "\\..*$", "")) %>%
    dplyr::filter(! is.na(log2FoldChange)) %>%
    left_join(dhs.reads.number.summarized.merged, by = c("DHS", "geno")) %>%
    dplyr::filter(scaled_total_per_lib >5) 

g <- ggplot(deseq.results.df_for.density, aes(y = log2FoldChange, fill = geno)) +
    geom_density() +
    theme_bw() +
    ## xlim(-1,  4) +
    scale_colour_manual(values = c("#3D405B","#F2CC8F")) +
    xlab("log2 Fold Change vs WT") +
    ggtitle("DHS log2FC upon exosome related factor KD")
ggsave(g, filename = file.path(out_dir, "DHS_log2FC_density_plot_more_than_5reads.pdf"), width = 8, height =4)

### establish list of DHS that are upregulated in each condition
deseq.results.df_up <- deseq.results.df %>%
    drop_na %>% 
    dplyr::filter(log2FoldChange > 1 & padj < 0.05  ) %>%
    dplyr::select(DHS, geno) %>%
    rbind.data.frame(data.frame(DHS = list.dhs.tested, geno = "WT")) %>%
    mutate(status = "upregulated") %>%
    remove_rownames

## make Venn diagram to see if the upregulated DHS are shared between the two conditions

venn_list <- split(deseq.results.df_up$DHS, deseq.results.df_up$geno)

venn_list$WT <- NULL 

g <- ggVennDiagram(venn_list) +
    scale_fill_gradient(low = "white", high = "red") +
    theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 2, unit = "cm")) +
    coord_cartesian(clip = "off")
ggsave(g, filename = file.path(out_dir, "DHS_upregulated_venn_diagram.pdf"), width = 6, height = 6)

list.dhs <- dhs.reads.number.summarized.merged %>%
    dplyr::select(DHS, geno, scaled_total_per_lib, strand) %>%
    mutate(DHS = paste(DHS, strand, sep = ".")) %>% 
    dplyr::filter(geno %in% c("WT",  "ZCCHC8",  "RBM7")) %>%
    dplyr::filter(scaled_total_per_lib > 1)%>%
    group_by( DHS) %>%
    add_count() %>%
    ungroup %>%
    dplyr::filter(n == 3) %>%
    .$DHS  %>% unique

deseq.results.df.xy <- deseq.results.df %>%
    dplyr::filter(DHS %in% list.dhs) %>%
    dplyr::select(DHS, log2FoldChange , geno) %>%
    reshape2::dcast(DHS ~ geno, value.var = "log2FoldChange", fill = 0)

custom_blues <- RColorBrewer::brewer.pal(9, "Blues")[3:9]

g <- ggplot(deseq.results.df.xy, aes(x=RBM7, y=ZCCHC8)) +
    geom_hex() +
    theme_bw() +
    geom_point(alpha = 0) +
    scale_fill_gradientn(colors = custom_blues, trans = "log10") +
    ggpubr::stat_cor(method = "pearson", label.x = -3, label.y = 7, size =  4 ) +
    xlab("log2 Fold Change RBM7 vs WT") +
    ylab("log2 Fold Change ZCCHC8 vs WT") +
    ggtitle("DHS log2FC upon exosome related factor KD") 
g <- ggMarginal(g, 
                type = "density",
                fill = "black")
ggsave(g, filename = file.path(out_dir, "DHS_log2FC_RBM7_vs_ZCCHC8_plot.pdf"), width = 6, height =4)

##### end of sensitivity plots We will also plot heatmaps later.
##############################

#### plot  length
dhs.up.read.length.pA.df <- dhs.reads.statistics.pA %>%
    dplyr::filter(geno != "ZC3H4") %>%
    mutate(geno = factor(geno, levels = c("WT",  "RBM7",  "ZCCHC8"))) %>% 
    dplyr::select(DHS, geno,  median.genomic.length, median.exonic.length, pA_status) %>% 
    left_join(deseq.results.df_up %>% separate(DHS, into=c("DHS","strand") , sep="\\."), by = c("DHS", "geno")) %>%
    dplyr::filter(!is.na(status))

g <- ggplot(dhs.up.read.length.pA.df, aes(x = median.exonic.length, colour = geno)) +
    geom_density() +
    xlim(0, 3000) +
    theme_bw() +
    facet_grid(rows = vars(pA_status)) +
    scale_colour_manual(values = c('#828282', "#3D405B","#F2CC8F")) +
    xlab("Median exonic read length") +
    ggtitle("DHS read genomic length distribution")
ggsave(g, filename = file.path(out_dir, "DHS_up_read_exonic_length_density_plot_pA.pdf"), width = 6, height =4)

#### plot exon number
dhs.up.read.length.pA.df <- dhs.reads.statistics.pA %>%
    mutate(geno = factor(geno, levels = c("WT",  "RBM7",  "ZCCHC8"))) %>%
    drop_na() %>% 
    left_join(deseq.results.df_up %>% separate(DHS, into=c("DHS","strand") , sep="\\."), by = c("DHS", "geno", "strand")) %>%
    drop_na() %>% 
    dplyr::select(DHS, geno, median.n.exon, pA_status) %>% 
    mutate(median.n.exon = round(median.n.exon + 0.1)) %>%  
    mutate(group_exon = case_when(
               median.n.exon == 1 ~ "1",
               median.n.exon == 2  ~ "2",
               median.n.exon >=3 & median.n.exon <=5 ~ "3-5",
               median.n.exon >=6 & median.n.exon <=10 ~ "6-10",
               median.n.exon >10 ~ ">10"
           )) %>%
    mutate(group_exon = factor(group_exon, levels = c("1", "2", "3-5", "6-10", ">10"))) %>%
    group_by(geno, pA_status) %>% 
    add_count() %>%
    ungroup() %>% 
    group_by(geno, pA_status, group_exon, n) %>%
    dplyr::summarize(count = n())  %>%
    mutate(fraction = count / n) %>%
    ungroup

g <- ggplot(dhs.up.read.length.pA.df, aes(x = group_exon, fill = geno, y = fraction)) +
    geom_histogram( stat = "identity", position = position_dodge()) +
    theme_bw() +
    facet_grid( rows = vars(pA_status)) +
    scale_fill_manual(values = c('#828282', "#3D405B","#F2CC8F")) +
    ggtitle("DHS read exon number fraction")
ggsave(g, filename = file.path(out_dir, "DHS_up_read_exon_number_fraction_histogram_pA.pdf"), width = 6, height =4)

### plot pA length
dhs.up.pA.length.df <- dhs.reads.statistics.pA %>%
    mutate(geno = factor(geno, levels = c("WT",  "RBM7",  "ZCCHC8"))) %>%
    drop_na() %>%
    dplyr::select(DHS, geno,  median.pA.length, pA_status)  %>%
    dplyr::filter(pA_status != "no_pA") %>%
    left_join(deseq.results.df_up %>% separate(DHS, into=c("DHS","strand") , sep="\\."), by = c("DHS", "geno")) %>%
    dplyr::filter(!is.na(status))

g <- ggplot(dhs.up.pA.length.df, aes(x = median.pA.length,  colour = geno)) +
    geom_density() +
    theme_bw() +
    scale_colour_manual(values = c('#828282', "#3D405B","#F2CC8F")) +
    ggtitle("DHS read pA length distribution") 
ggsave(g, filename = file.path(out_dir, "DHS_read_pA_length_density_plot.pdf"), width = 6, height =4)

### plot pA fraction (all expressed PROMPTs, not only upregulated)

dhs.pA.fraction.df <- dhs.reads.number.summarized.pA %>%
    dplyr::filter(geno %in% c("WT",  "RBM7",  "ZCCHC8")) %>%
    dplyr::select(DHS, geno ,  scaled_total_per_lib, pA_status) %>%
    reshape2::dcast(DHS + geno  ~ pA_status , value.var = "scaled_total_per_lib", fill = 0) %>%
    mutate(fraction = pA / (pA + no_pA)) %>%
    left_join(dhs.reads.number.summarized.merged, by = c("DHS", "geno")) %>%
    dplyr::filter(scaled_total_per_lib > 5) %>% 
    mutate(geno = factor(geno, levels = c("WT",  "RBM7",  "ZCCHC8"))) 

g <- ggplot(dhs.pA.fraction.df, aes(x = fraction, colour = geno, fill = geno,
                                    y = after_stat(count))) +
    geom_histogram(binwidth = 0.05) +
    theme_bw() +
    facet_grid(row = vars(geno)) +
    scale_colour_manual(values = c('#828282', "#3D405B","#F2CC8F")) +
    scale_fill_manual(values = c('#828282', "#3D405B","#F2CC8F")) +
    xlab("pA read fraction") +
    ggtitle("DHS read pA fraction distribution")
ggsave(g, filename = file.path(out_dir, "DHS_read_pA_fraction_histogram_same_scale.pdf"), width = 6, height =4)

#### same, but only up-regulated (not for paper)

dhs.up.pA.fraction.df <- dhs.reads.number.summarized.pA %>%
    dplyr::filter(geno %in% c("WT",  "RBM7",  "ZCCHC8")) %>%
    dplyr::select(DHS, geno ,  scaled_total_per_lib, pA_status) %>%
    reshape2::dcast(DHS + geno  ~ pA_status , value.var = "scaled_total_per_lib", fill = 0) %>%
    mutate(fraction = pA / (pA + no_pA)) %>%
    left_join(dhs.reads.number.summarized.merged, by = c("DHS", "geno")) %>%
    dplyr::filter(scaled_total_per_lib > 5) %>% 
    left_join(deseq.results.df_up %>% separate(DHS, into=c("DHS","strand") , sep="\\."), by = c("DHS", "geno")) %>%
    dplyr::filter(!is.na(status))  %>%
    mutate(geno = factor(geno, levels = c("WT",  "RBM7",  "ZCCHC8")))  

g <- ggplot(dhs.up.pA.fraction.df, aes(x = fraction, colour = geno, fill = geno,
                                    y = after_stat(count))) +
    geom_histogram(binwidth = 0.05) +
    theme_bw() +
    facet_grid(row = vars(geno)) +
    scale_colour_manual(values = c('#828282', "#3D405B","#F2CC8F")) +
    scale_fill_manual(values = c('#828282', "#3D405B","#F2CC8F")) +
    xlab("pA read fraction") +
    ggtitle("DHS read pA fraction distribution")
ggsave(g, filename = file.path(out_dir, "DHS_up_read_pA_fraction_histogram_same_scale.pdf"), width = 6, height =4)

#### make violin plots pA fraction
##################################

dhs.pA.fraction.df <- dhs.reads.number.summarized.pA %>%
    dplyr::filter(geno %in% c("WT",   "RBM7",  "ZCCHC8")) %>%
    dplyr::select(DHS, geno, strand,  scaled_total_per_lib, pA_status)   %>%
    reshape2::dcast(DHS + geno + strand ~ pA_status, value.var = "scaled_total_per_lib", fill = 0) %>%
    mutate(fraction = pA / (pA + no_pA)) %>%
    left_join(dhs.reads.number.summarized.merged, by = c("DHS", "geno", "strand")) %>%
    dplyr::filter(scaled_total_per_lib > 5) %>% 
    mutate(geno = factor(geno, levels = c("WT", "RBM7",  "ZCCHC8")))

dhs.read.length.merged.df <- dhs.reads.statistics.merged %>%
    dplyr::filter(geno %in% c("WT",   "RBM7",  "ZCCHC8")) %>%
    mutate(geno = factor(geno, levels = c("WT",   "RBM7",  "ZCCHC8"))) %>%
    dplyr::select(DHS, geno, strand, median.exonic.length)

dhs.list.expressed <- dhs.reads.number.summarized.merged  %>%
    dplyr::filter(scaled_total_per_lib > 5) %>%
    dplyr::filter(geno %in% c("WT",   "RBM7",  "ZCCHC8"))
write.table(dhs.list.expressed %>% dplyr::select(DHS, strand) %>% distinct,
           file = file.path(out_dir, "DHS_list_expressed_reads_over_5.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

logFC.readLength.df <- left_join(dhs.read.length.merged.df,
                                 deseq.results.df %>% separate(DHS, into=c("DHS", "strand"), sep = "\\."),
                                 by = c("DHS", "strand", "geno")) %>%
    left_join(dhs.pA.fraction.df %>% dplyr::select(DHS, strand, geno, fraction),
              by = c("DHS", "strand", "geno")) %>%
    dplyr::select(-padj, -metric)  %>%
    mutate(fraction = ifelse(is.na(fraction), 0, fraction)) %>%
    drop_na %>%
    dplyr::rename(fractionpA = fraction) %>%
    left_join(dhs.list.expressed, by=c("DHS","geno","strand")) %>%
    drop_na() %>%
    mutate(geno = factor(geno, levels = c(  "RBM7",   "ZCCHC8")))

## make intervals

cI <- classIntervals(logFC.readLength.df$log2FoldChange, n = 6, style = "kmeans")
kmeans_breaks <- round(cI$brks, 2)
kmeans_breaks[1] <- kmeans_breaks[1] - 0.01
kmeans_breaks[length(kmeans_breaks)] <- kmeans_breaks[length(kmeans_breaks)] + 0.01
left_bounds <- kmeans_breaks[-length(kmeans_breaks)]
right_bounds <- kmeans_breaks[-1]
interval_strings <- sprintf("[%s, %s)", left_bounds, right_bounds)
num_intervals <- length(interval_strings)
interval_strings[num_intervals] <- sprintf("[%s, %s]",  left_bounds[num_intervals],  right_bounds[num_intervals])
indices <- findInterval(logFC.readLength.df$log2FoldChange, kmeans_breaks, rightmost.closed = TRUE)
interval_assignments <- factor(interval_strings[indices], levels = interval_strings)
saveRDS(cI, file.path(out_dir, "CI_prompts"))

logFC.readLength.df <- logFC.readLength.df %>%
    mutate(id3 = interval_assignments)

## prepare the plot
counts <- logFC.readLength.df %>%
    group_by(geno, id3) %>%
    summarise(count = n()) %>%
    mutate(count  = paste0("n = ", count)) %>%
    mutate(fractionpA= 0)

ano.list <- list()
for (geno.elt in unique(logFC.readLength.df$geno))
{
    ano.list.df <- logFC.readLength.df %>%
        dplyr::filter(geno.elt == geno)
    ano.list[[geno.elt]] <- lm(fractionpA ~ id3, data = ano.list.df) %>% anova
}

ano.df <- ano.list %>% sapply(getElement, "Pr(>F)") %>% .[1,] %>% data.frame(geno = names(.), pval = .) %>%
    mutate(cat = paste("ANOVA pval =", formatC(pval, format = "e", digits = 2))) %>%
    remove_rownames() %>%
    dplyr::select(-pval)

p1 <- ggplot(logFC.readLength.df,  aes(y = fractionpA, x = id3, fill = geno)) +
    geom_violin( color = "black", size = 0.2) +
    theme_bw() +
    theme(# axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x  = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.placement = "outside") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
    stat_summary(fun =mean, geom="point", shape=1, size=2, color="black") +
    stat_summary(fun = median, geom="point", shape=2, size=2, color="black") +
    facet_grid(row = vars(geno),  switch = "x") +
    scale_fill_manual(values = c( "#3D405B","#F2CC8F")) +
    ylab("fraction pA reads") +
    ylim(c(0,1)) +
    xlab("log2FC (in k-mean intervals)") +
    geom_text(aes(label = cat), data = logFC.readLength.df %>% left_join(ano.df) %>% distinct(cat, .keep_all = TRUE),
              size=3, y = 0.95, x = "[1.66, 2.6)")
ggsave(p1, filename = file.path(out_dir, "DHS_fraction_pA_vs_logFC_violin_plot_2.pdf"), width = 6, height = 4)

########################################
### make heatmaps sensityivity vs length.

deseq.results.df.pheatmap <- deseq.results.df %>%
    dplyr::filter(DHS %in% list.dhs) %>%
    dplyr::select(DHS, log2FoldChange , geno) %>%
    mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange) %>%
               ifelse(. > 6 , 6, .) %>% ifelse(. < -3, -3, .)) %>%
    reshape2::dcast(DHS ~ geno, value.var = "log2FoldChange", fill = 0) %>%
    column_to_rownames("DHS")

dist.mat <- dist(deseq.results.df.pheatmap)
h <- hclust(dist.mat)

DHS.order <- h$labels[h$order]

col.pal <-  c(colorRampPalette(rev(brewer.pal(n = 9, name = "PuRd"))[c(4,5,6,7,9)])(33), "white",
              colorRampPalette(rev(brewer.pal(n = 9, name = "Blues"))[3:9])(66) %>% rev)

## col.pal <- c(colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu"))[1:4])(33),
##              colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu"))[4:7])(66))

p <- pheatmap(deseq.results.df.pheatmap[DHS.order, c( "RBM7", "ZCCHC8")], 
              cluster_rows = FALSE, 
              cluster_cols = FALSE,
              show_rownames = FALSE,
              fontsize_col = 10,
              breaks = seq(-3,6, length.out = 100),
              color = col.pal,
              width = 6, height = 10)
ggsave(file.path(out_dir,  "sensitivity_heatmap.pdf"), p, width=4, height=10)

dhs.read.length.merged.df <- dhs.reads.statistics.merged %>%
    dplyr::filter(geno %in% c("WT","RBM7","ZCCHC8")) %>%
    mutate(geno = factor(geno, levels = c("WT", "RBM7", "ZCCHC8"))) %>%
    dplyr::select(DHS, geno, strand, median.exonic.length)

dhs.read.length.merged.pheatmap <- dhs.read.length.merged.df %>%
    mutate(DHS = paste(DHS, strand, sep = ".")) %>%
    mutate(median.exonic.length = ifelse(median.exonic.length > 1000, 1000, median.exonic.length)) %>%
    reshape2::dcast(DHS  ~ geno, value.var = "median.exonic.length", fill = 0)  %>%
    column_to_rownames("DHS")

col.pal.2 <-  c(colorRampPalette(rev(brewer.pal(n = 9, name = "PuRd"))[c(4,5,6,7,9)])(50), 
              colorRampPalette(rev(brewer.pal(n = 9, name = "Blues"))[3:9])(50) %>% rev) %>% rev
p <- pheatmap(dhs.read.length.merged.pheatmap[DHS.order, c("WT", "RBM7", "ZCCHC8")], 
              cluster_rows = FALSE, 
              cluster_cols = FALSE,
              show_rownames = FALSE,
              fontsize_col = 10,
              width = 6, height = 10,
              color = col.pal.2)
ggsave(file.path(out_dir,  "median_exonic_length_heatmap.pdf"), p, width=4, height=10)


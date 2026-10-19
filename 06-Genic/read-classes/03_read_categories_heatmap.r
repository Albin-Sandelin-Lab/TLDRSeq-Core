## Frequency of the read categories (3' end) per gene-associated DHS: heatmap in WT, and
## fraction of over-elongated reads.
##
## Usage (from 06-Genic/read-classes/): Rscript 03_read_categories_heatmap.r
##
## Input : GENIC_FREQ_RDS (config/paths.R; ../proportions/02_build_freq_table.R): frequency of each
##         category per DHS and genotype (TLDRSeq_sup)
## Output (in $GENIC_DIR/TLDRSeq_sup/read_categories/):
##   pheatmap_scaled_freq_WT.pdf               WT frequencies, DHSs ordered by clustering (Ward)
##   over_elongation_fraction.tsv              fraction of over-elongated reads in WT
##
## Frequencies: read counts per DHS, replicate and category (pA+ / pA- separated) divided
## by the chrM scale factor of the replicate and summed over the replicates of each genotype;
## frequency = category count / all reads of the DHS in the genotype.

library(tidyverse)
library(pheatmap)

source("../../config/paths.R")

dataset <- "TLDRSeq_sup"
out_dir <- file.path(GENIC_DIR, dataset, "read_categories")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genic.freq <- readRDS(GENIC_FREQ_RDS)


#### Heatmap of the category frequencies in WT

tab.wide <- genic.freq %>%
    dplyr::filter(geno %in% c("WT")) %>%
    dplyr::filter(feature_id %in% c("typical_TTS_pAp", "typical_TTS_pAn",
                                    "gradual_degradation_pAp", "gradual_degradation_pAn",
                                    "first_exon_pAp",  "first_exon_pAn",
                                    "premature_termination_pAp", "premature_termination_pAn")) %>%
    mutate(feature_id = factor(feature_id, c("typical_TTS_pAp", "typical_TTS_pAn",
                                             "gradual_degradation_pAp", "gradual_degradation_pAn",
                                             "premature_termination_pAp", "premature_termination_pAn",
                                             "first_exon_pAp",  "first_exon_pAn"))) %>%
    reshape2::dcast(DHS ~ feature_id + geno, value.var = "scaled_freq", fill = 0) %>%
    column_to_rownames("DHS")

d <- dist(tab.wide)
order.d <- hclust(d, "ward.D")$order

col_fun_prop <- colorRampPalette(c("white", "#028090"))(100)
breaks_prop <- seq(0, 1, length.out = 101)

pheatmap(tab.wide[order.d, ] %>% t,
         color = col_fun_prop, breaks = breaks_prop,
         show_rownames = TRUE, show_colnames = FALSE,
         cluster_rows = FALSE, cluster_cols = FALSE,
         filename = file.path(out_dir, "pheatmap_scaled_freq_WT.pdf"),
         width = 8, height = 4)


#### Fraction of over-elongated reads in WT

genic.freq.WT <- genic.freq %>%
    dplyr::filter(geno %in% c("WT"))

## all DHSs pooled: over-elongated reads / all reads except "others"
n.over_elong.c <- genic.freq.WT %>%
    dplyr::filter(feature_id %in% c("over_elongation_pAp", "over_elongation_pAn")) %>%
    .$feature_count %>%
    sum

n.not.over_elong.c <- genic.freq.WT %>%
    dplyr::filter(!feature_id %in% c("others")) %>%
    .$feature_count %>%
    sum

frac.pooled <- n.over_elong.c / n.not.over_elong.c

## per DHS, then averaged over DHSs
frac.mean.per.DHS <- genic.freq.WT %>%
    dplyr::filter(feature_id != "others") %>%
    mutate(annot = ifelse(feature_id %in% c("over_elongation_pAp", "over_elongation_pAn"), "num", "rest")) %>%
    group_by(DHS, annot) %>%
    summarise(feature_count = sum(feature_count), .groups = "drop") %>%
    group_by(DHS) %>%
    mutate(total = sum(feature_count)) %>%
    dplyr::filter(annot == "num") %>%
    dplyr::filter(total != 0) %>%
    mutate(prop = feature_count / total) %>%
    .$prop %>%
    mean()

df.over_elong <- data.frame(measure = c("all DHSs pooled", "mean over DHSs"),
                            fraction_over_elongation = c(frac.pooled, frac.mean.per.DHS))
print(df.over_elong)
write.table(df.over_elong, file.path(out_dir, "over_elongation_fraction.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

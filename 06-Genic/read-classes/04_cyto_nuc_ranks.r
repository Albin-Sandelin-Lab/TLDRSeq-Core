## Cytoplasmic / nuclear localization of the read categories (TLDRSeq_nuclear_sup): DHS x
## category ranked by their cytoplasmic / nuclear ratio.
##
## Usage (from 06-Genic/read-classes/): Rscript 04_cyto_nuc_ranks.r
##
## Input : $GENIC_DIR/TLDRSeq_nuclear_sup/attribute_reads/df_unified_results_classified.rds (01)
##         $GENIC_DIR/TLDRSeq_nuclear_sup/library_size_factors/library_size_factors.tsv     (02)
## Output (in $GENIC_DIR/TLDRSeq_nuclear_sup/cyto_nuc_ranks/):
##   heatmap_rank_cyto_nuc_by_3prime_category.png      rank of each gene in 4 categories
##                                                     (genes with reads in all 4)
##   heatmap_rank_cyto_nuc_by_3prime_category_all.png  same, all categories (pA+ and pA-)
##   tab_wide_rank_long.tsv       ranks, 4 categories, DHSs with reads in all 4, with deciles
##   tab_wide_rank_long.full.tsv  ranks, 4 categories (pA+ and pA-), with deciles per category
##   tab_wide_rank_first_exon.tsv ranks of first exonic pA- reads (rank among all categories), deciles
##
## Ratio: (mean cytoplasmic CPM + 1) / (mean nuclear CPM + 1) over the replicates;
## rank 1 = most cytoplasmic. CPM = read count / library depth * 1e6.
## No scaling between fractions (no spike-in; no mitochondrial RNA in the nuclear fraction),
## hence ranks. The total RNA libraries ("full") are not used.

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

source("../../config/paths.R")

dataset <- "TLDRSeq_nuclear_sup"
in_dir  <- file.path(GENIC_DIR, dataset)
out_dir <- file.path(GENIC_DIR, dataset, "cyto_nuc_ranks")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

object <- readRDS(file.path(in_dir, "attribute_reads", "df_unified_results_classified.rds"))

## sample (e.g. nucl_1, Cyto_2), depth
norm.table <- read.table(file.path(in_dir, "library_size_factors", "library_size_factors.tsv"), header = TRUE, sep = "\t")


#### CPM per DHS, replicate and category

## upstream, over-elongated and other reads are not used (also not for the ranks)
df_DHS_prop_3prime_end_curated.base <- object %>%
    dplyr::filter(!label_3prime %in% c("Others", "upstream", "over_elongation")) %>%
    group_by(sample) %>%
    dplyr::filter(strand.x == strand.y) %>%
    dplyr::filter(exonCount != 1, !grepl("full", sample)) %>%
    group_by(gene_name, DHS, sample, strand.x, label_3prime_pA, polyA_status, label_3prime) %>%
    summarise(
        n = n()
    ) %>%
    left_join(norm.table %>% dplyr::select(sample, depth)) %>%
    mutate(CPM = n / depth * 1e6)


#### Cytoplasmic / nuclear ratio

df_ratios <- df_DHS_prop_3prime_end_curated.base %>%
    ungroup %>%
    separate(sample, into = c("fraction", "replicate"), sep = "_") %>%
    group_by(fraction, gene_name, DHS, label_3prime_pA, polyA_status, label_3prime) %>%
    dplyr::summarize(n = mean(CPM)) %>%
    ungroup() %>%
    reshape2::dcast(gene_name + DHS + label_3prime_pA + polyA_status + label_3prime ~ fraction, value.var = "n", fill = 0) %>%
    mutate(
        cyto_nuc_ratio = (Cyto + 1) / (nucl + 1),
        log2_ratio = log2(cyto_nuc_ratio)
    ) %>%
    arrange(desc(log2_ratio)) %>%
    mutate(overall_rank = row_number()) %>%
    dplyr::select(gene_name, DHS, label_3prime_pA, log2_ratio, overall_rank, label_3prime, polyA_status, nucl, Cyto, cyto_nuc_ratio) %>%
    mutate(log10_rank = log10(overall_rank + 1))


#### Heatmap: rank of each gene in the 4 categories

tab_wide.rank <- df_ratios %>%
    ungroup() %>%
    mutate(value = cyto_nuc_ratio) %>%
    dplyr::filter(label_3prime_pA %in%
                  c("premature_termination_pA-", "first_exon_pA-", "typical_TTS_pA+", "gradual_degradation_pA-")) %>%
    arrange(desc(value)) %>%
    mutate(overall_rank = row_number()) %>%
    mutate(overall_rank.norm = overall_rank / nrow(.)) %>%
    mutate(log_overall_rank = log10(overall_rank)) %>%
    reshape2::dcast(gene_name ~ label_3prime_pA, value.var = "overall_rank", fill = NA) %>%
    column_to_rownames("gene_name") %>%
    setNames(colnames(.) %>% str_replace_all("_", " ")) %>%
    drop_na %>%
    setNames(c("First exonic (pA-)", "Intermediate (pA-)", "First intronic (pA-)", "Full-length (pA+)"))

d <- dist(cbind(tab_wide.rank,
                tab_wide.rank[, 1],
                tab_wide.rank[, 3]
                ))
order.d <- hclust(d)$order %>% rev

custom_blues <- RColorBrewer::brewer.pal(9, "Blues")[1:9]
p <- pheatmap(tab_wide.rank[order.d, 4:1] %>% t,
              color = colorRampPalette(custom_blues)(100),
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              show_colnames = FALSE,
              show_rownames = TRUE,
              fontsize = 8,
              filename = file.path(out_dir, "heatmap_rank_cyto_nuc_by_3prime_category.png"),
              width = 5,
              height = 2
              )


#### Heatmap: all categories, pA+ and pA- (NA = no reads)

tab_wide.rank.all <- df_ratios %>%
    ungroup() %>%
    mutate(value = cyto_nuc_ratio) %>%
    arrange(desc(value)) %>%
    mutate(overall_rank = row_number()) %>%
    mutate(overall_rank.norm = overall_rank / nrow(.)) %>%
    mutate(log_overall_rank = log10(overall_rank)) %>%
    reshape2::dcast(gene_name ~ label_3prime_pA, value.var = "overall_rank", fill = 0) %>%
    column_to_rownames("gene_name") %>%
    setNames(colnames(.) %>% str_replace_all("_", " ")) %>%
    dplyr::rename(
               "First exonic (pA-)" = "first exon pA-",
               "First exonic (pA+)" = "first exon pA+",
               "Intermediate (pA-)" = "gradual degradation pA-",
               "Intermediate (pA+)" = "gradual degradation pA+",
               "First intronic (pA-)" = "premature termination pA-",
               "First intronic (pA+)" = "premature termination pA+",
               "Full-length (pA+)" = "typical TTS pA+",
               "Full-length (pA-)" = "typical TTS pA-"
           )

d <- dist(cbind(tab_wide.rank.all,
                tab_wide.rank.all[, 5]
                ))

order.d <- hclust(d)$order

tab_wide.rank.all[tab_wide.rank.all == 0] <- NA

p <- pheatmap(tab_wide.rank.all[order.d, ] %>% t,
              color = colorRampPalette(custom_blues)(100),
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              show_colnames = FALSE,
              fontsize = 8,
              filename = file.path(out_dir, "heatmap_rank_cyto_nuc_by_3prime_category_all.png"),
              width = 15,
              height = 8
              )


#### Rank tables

## 4 categories, DHSs with reads in all 4; deciles of rank per category
tab_wide.rank.long <- df_ratios %>%
    ungroup() %>%
    dplyr::filter(label_3prime_pA %in%
                  c("premature_termination_pA-", "first_exon_pA-", "typical_TTS_pA+", "gradual_degradation_pA-")) %>%
    mutate(value = cyto_nuc_ratio) %>%
    arrange(desc(value)) %>%
    mutate(overall_rank = row_number()) %>%
    mutate(overall_rank.norm = overall_rank / nrow(.)) %>%
    mutate(log_overall_rank = log10(overall_rank)) %>%
    reshape2::dcast(DHS ~ label_3prime_pA, value.var = "overall_rank", fill = NA) %>%
    column_to_rownames("DHS") %>%
    drop_na %>%
    rownames_to_column("DHS") %>%
    reshape2::melt(id.vars = "DHS", variable.name = "label_3prime_pA", value.name = "overall_rank") %>%
    drop_na() %>%
    group_by(label_3prime_pA) %>%
    mutate(nqtile.per.label = ntile(overall_rank, 10))

write.table(tab_wide.rank.long, file.path(out_dir, "tab_wide_rank_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## 4 categories, pA+ and pA-, all DHSs; deciles of rank per category (used by 05)
tab_wide.rank.long.full <- df_ratios %>%
    dplyr::filter(label_3prime %in%
                  c("premature_termination", "first_exon", "typical_TTS", "gradual_degradation")) %>%
    ungroup() %>%
    mutate(value = cyto_nuc_ratio) %>%
    arrange(desc(value)) %>%
    mutate(overall_rank = row_number()) %>%
    group_by(label_3prime_pA) %>%
    mutate(nqtile.per.label = ntile(overall_rank, 10))

write.table(tab_wide.rank.long.full, file.path(out_dir, "tab_wide_rank_long.full.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## first exonic pA- reads: rank among all DHS x categories, with deciles (used by the genic AT models)
tab_wide.rank.first_exon.df <- df_ratios %>%
    ungroup() %>%
    dplyr::filter(label_3prime_pA %in% c("first_exon_pA-")) %>%
    mutate(value = cyto_nuc_ratio) %>%
    arrange(desc(value)) %>%
    mutate(nqtile = ntile(overall_rank, 10)) %>%
    dplyr::select(DHS, label_3prime_pA, overall_rank, nqtile)

write.table(tab_wide.rank.first_exon.df, file.path(out_dir, "tab_wide_rank_first_exon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

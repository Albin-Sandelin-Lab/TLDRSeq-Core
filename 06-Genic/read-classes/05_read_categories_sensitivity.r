## Sensitivity of the read categories to the exosome cofactors: pA+ fraction, cytoplasmic /
## nuclear rank and DESeq2 log2FC per category, for NEXT (RBM7, ZCCHC8) and PAXT (MTR4, Z1-FL)
## samples.
##
## Usage (from 06-Genic/read-classes/): Rscript 05_read_categories_sensitivity.r
##
## Input : TTS_PROP_CATEGORIES_RDS (config/paths.R; label_reads_TLDRSeq_sup.R): TLDRSeq_sup reads
##           per DHS, replicate and 3' end label combination
##         $GENIC_DIR/TLDRSeq_Z1_MTR4_sup/attribute_reads/df_unified_results_classified.rds (01)
##         $GENIC_DIR/TLDRSeq_Z1_MTR4_sup/library_size_factors/library_size_factors.tsv (02)
##         $GENIC_DIR/TLDRSeq_nuclear_sup/cyto_nuc_ranks/tab_wide_rank_long.full.tsv     (04)
##         chrM scale factors of TLDRSeq_sup: SCALE_FACTORS_CSV, CHRM_DEPTH (config/paths.R)
## Output (in $GENIC_DIR/read_categories_sensitivity/):
##   fraction_pA_plus_by_3prime_category.pdf          fraction of pA+ reads per DHS and category
##                                                     (WT, RBM7, ZCCHC8)
##   rank_of_DHS_by_3prime_category.pdf               cytoplasmic / nuclear rank of the DHS x
##                                                     category detected in WT, RBM7, ZCCHC8
##   deseq_log2FC_by_3prime_category_merged.pdf       DESeq2 log2FC (RBM7, ZCCHC8 vs WT) per category
##   violin_log2FoldChange_NEXT.pdf                   log2FC of first exonic pA- reads of the most
##                                                     cytoplasmic vs most nuclear DHSs
##   deseq_log2FC_by_3prime_category_Z1_MTR4.pdf      DESeq2 log2FC (MTR4, Z1-FL vs WT)
##
## Categories: Full-length = typical_TTS, Intermediate = gradual_degradation,
## First intronic = premature_termination, First exonic = first_exon.

library(tidyverse)
library(DESeq2)

source("../../config/paths.R")
source("../../lib/tldr_reads.R")
source("../../lib/deseq.R")
source("../../04-Normalization/scale_factors.R")

out_dir <- file.path(GENIC_DIR, "read_categories_sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

main_categories <- c("premature_termination", "first_exon", "typical_TTS", "gradual_degradation")

## Read counts per DHS, replicate and category (pA+ / pA- separated)
## TLDRSeq_sup: categories from the label combinations (label_reads_TLDRSeq_sup.R)
count_reads_per_category_sup <- function() {
    readRDS(TTS_PROP_CATEGORIES_RDS) %>%
        mutate(label_3prime = case_when(
                   label_tts_falling_region == "upstream" ~ "upstream",
                   label_premature_termination == "singleExon" ~ "single_exon",
                   label_tts == "before_TTS" & label_premature ~ "premature_termination",
                   label_tts == "before_TTS" & label_premature_termination == "firstExon" & unspliced ~ "first_exon",
                   label_tts == "before_TTS" & label_cdsts == "before_CDS" &
                   label_premature_termination == "others" ~ "gradual_degradation",
                   label_tts == "before_TTS" & (label_cdsts == "after_CDS" | label_cdsts == "at_CDS") ~ "UTR3",
                   label_tts == "typical_TTS" ~ "typical_TTS",
                   label_tts == "overlong_TTS" ~ "over_elongation",
                   TRUE ~ "Others"
               )) %>%
        dplyr::filter(label_3prime %in% main_categories) %>%
        mutate(label_3prime_pA = paste0(label_3prime, "_", ifelse(polyA_tail, "pA+", "pA-")),
               pA_status = ifelse(polyA_tail, "pA", "no_pA")) %>%
        group_by(DHS, gene_name, sample, label_3prime, label_3prime_pA, pA_status) %>%
        summarise(n = sum(n), .groups = "drop")
}

## TLDRSeq_Z1_MTR4_sup: categories of 01_classify_reads_per_DHS.R
count_reads_per_category <- function(dataset) {
    readRDS(file.path(GENIC_DIR, dataset, "attribute_reads", "df_unified_results_classified.rds")) %>%
        mutate(pA_status = ifelse(polyA_status == "pA+", "pA", "no_pA")) %>%
        dplyr::count(DHS, gene_name, sample, label_3prime, label_3prime_pA, pA_status)
}


#### Load data

## chrM scale factors (TLDRSeq_sup)
df_scale_factor <- read.csv(SCALE_FACTORS_CSV) %>%
    mutate(Rep = sample %>% str_replace("scale_", "")) %>%
    dplyr::select(scale_factor, sample = Rep)

## read counts per category (TLDRSeq_sup)
df_counts <- count_reads_per_category_sup()

## cytoplasmic vs nuclear ranks (TLDRSeq_nuclear_sup)
dhs.cyto.vs.nuc <- read.table(file.path(GENIC_DIR, "TLDRSeq_nuclear_sup", "cyto_nuc_ranks", "tab_wide_rank_long.full.tsv"), sep = "\t", header = TRUE) %>%
    dplyr::select(DHS, overall_rank, label_3prime_pA) %>%
    mutate(normalized_rank = overall_rank / max(overall_rank))

tab.categories <- read.table(file.path(GENIC_DIR, "TLDRSeq_nuclear_sup", "cyto_nuc_ranks", "tab_wide_rank_long.full.tsv"), sep = "\t", header = TRUE) %>%
    dplyr::select(DHS, nqtile.per.label, label_3prime_pA, overall_rank)


#### Scaled read counts per DHS, genotype and category

tab.classified <- df_counts %>%
    dplyr::filter(label_3prime %in% main_categories) %>%
    left_join(df_scale_factor, by = "sample") %>%
    mutate(n_scaled = n / scale_factor) %>%
    separate(sample, into = c("garb", "geno", "garb2"), sep = "_") %>%
    dplyr::select(-garb, -garb2)


#### Fraction of pA+ reads per category

## mean scaled counts per replicate (6 WT replicates, 3 per depletion); DHSs with > 5 reads
tab.classified.to.plot <- tab.classified %>%
    dplyr::filter(geno %in% c("WT", "RBM7", "ZCCHC8")) %>%
    group_by(DHS, gene_name, label_3prime, label_3prime_pA, pA_status, geno) %>%
    summarise(n_scaled = sum(n_scaled)) %>%
    mutate(n_scaled = ifelse(geno == "WT", n_scaled / 6, n_scaled / 3)) %>%
    group_by(DHS, gene_name, geno) %>%
    mutate(n_scaled_tot_pergene = sum(n_scaled)) %>%
    dplyr::filter(n_scaled_tot_pergene > 5) %>%
    reshape2::dcast(gene_name + DHS + geno + label_3prime ~ pA_status, value.var = "n_scaled", fill = 0) %>%
    mutate(fraction = pA / (pA + no_pA)) %>%
    mutate(geno = factor(geno, levels = c("WT", "RBM7", "ZCCHC8"))) %>%
    mutate(label_3prime = case_when(
               label_3prime == "first_exon" ~ "First exonic",
               label_3prime == "gradual_degradation" ~ "Intermediate",
               label_3prime == "premature_termination" ~ "First intronic",
               label_3prime == "typical_TTS" ~ "Full-length")
           ) %>%
    mutate(label_3prime = factor(label_3prime,
                                 levels = c("Full-length", "Intermediate", "First intronic", "First exonic")))

p <- ggplot(tab.classified.to.plot,
            aes(x = label_3prime, y = fraction, fill = geno)) +
    geom_violin(scale = "width", colour = "transparent") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "3' end category", y = "Fraction of pA+ reads") +
    scale_fill_manual(values = c('#7b7c79', "#b1cdc3", "#886e97"))
ggsave(p, filename = file.path(out_dir, "fraction_pA_plus_by_3prime_category.pdf"), width = 10, height = 4)


#### Cytoplasmic / nuclear rank of the DHS x category detected in each genotype

tab.classified.to.plot.2 <- tab.classified %>%
    dplyr::filter(geno %in% c("WT", "RBM7", "ZCCHC8")) %>%
    group_by(DHS, gene_name, label_3prime, label_3prime_pA, pA_status, geno) %>%
    summarise(n_scaled = sum(n_scaled)) %>%
    mutate(n_scaled = ifelse(geno == "WT", n_scaled / 6, n_scaled / 3)) %>%
    group_by(DHS, gene_name, geno) %>%
    mutate(n_scaled_tot_pergene = sum(n_scaled)) %>%
    dplyr::filter(n_scaled_tot_pergene > 5) %>%
    dplyr::select(-n_scaled_tot_pergene) %>%
    left_join(dhs.cyto.vs.nuc, by = c("DHS", "label_3prime_pA"), multiple = "all") %>%
    dplyr::select(-label_3prime_pA) %>%
    dplyr::filter(!is.na(overall_rank)) %>%
    mutate(geno = factor(geno, levels = c("WT", "RBM7", "ZCCHC8"))) %>%
    mutate(label_3prime = case_when(
               label_3prime == "first_exon" ~ "First exonic",
               label_3prime == "gradual_degradation" ~ "Intermediate",
               label_3prime == "premature_termination" ~ "First intronic",
               label_3prime == "typical_TTS" ~ "Full-length")
           ) %>%
    mutate(label_3prime = factor(label_3prime,
                                 levels = c("Full-length", "Intermediate", "First intronic", "First exonic")))

g <- ggplot(tab.classified.to.plot.2, aes(x = label_3prime, y = normalized_rank, fill = geno)) +
    geom_violin(scale = "width", colour = "transparent") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "3' end category", y = " Normalized Rank of DHS * label (cyto vs nuc)") +
    scale_fill_manual(values = c('#7b7c79', "#b1cdc3", "#886e97"))
ggsave(g, filename = file.path(out_dir, "rank_of_DHS_by_3prime_category.pdf"), width = 10, height = 4)


#### DESeq2: RBM7, ZCCHC8 vs WT, per DHS x category

tab.deseq <- df_counts %>%
    dplyr::filter(label_3prime %in% main_categories) %>%
    dplyr::select(DHS, gene_name, sample, label_3prime_pA, n) %>%
    group_by(DHS, gene_name, sample, label_3prime_pA) %>%
    summarise(n = sum(n)) %>%
    mutate(DHS.label = paste0(DHS, ".", label_3prime_pA)) %>%
    dplyr::filter(!grepl("ZC3H4|ZC3H18|ARS2", sample)) %>%
    reshape2::dcast(DHS.label ~ sample, value.var = "n", fill = 0) %>%
    column_to_rownames("DHS.label")

df.meta <- colnames(tab.deseq) %>%
    data.frame(sample = .) %>%
    mutate(geno = sample %>% str_replace("mES_", "") %>% str_replace("_rep.", "")) %>%
    column_to_rownames("sample") %>%
    mutate(geno = factor(geno, c("WT", "RBM7", "ZCCHC8")))

deseq.results <- run_deseq(tab.deseq, df.meta, candidate_geno = c("RBM7", "ZCCHC8"))

deseq.results.df <- deseq.results$df_lfc_shrunk %>%
    separate(DHS, into = c("garb", "label_3prime_pA"), sep = "\\.") %>%
    dplyr::select(-garb) %>%
    mutate(label_3prime = label_3prime_pA %>% str_replace("_pA\\+", "") %>% str_replace("_pA-", "")) %>%
    mutate(pA_status = ifelse(grepl("_pA\\+", label_3prime_pA), "pA", "no_pA")) %>%
    dplyr::select(-label_3prime_pA) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    dplyr::select(label_3prime, pA_status, log2FoldChange, geno) %>%
    mutate(label_3prime = case_when(
               label_3prime == "first_exon" ~ "First exonic",
               label_3prime == "gradual_degradation" ~ "Intermediate",
               label_3prime == "premature_termination" ~ "First intronic",
               label_3prime == "typical_TTS" ~ "Full-length")
           ) %>%
    mutate(label_3prime = factor(label_3prime,
                                 levels = c("Full-length", "Intermediate", "First intronic", "First exonic")))

g <- ggplot(deseq.results.df, aes(x = label_3prime, y = log2FoldChange, fill = geno)) +
    geom_violin(scale = "width", colour = "transparent") +
    theme_bw() +
    coord_cartesian(ylim = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "3' end category", y = "log2 Fold Change (vs WT)") +
    scale_fill_manual(values = c("#b1cdc3", "#886e97"))
ggsave(g, filename = file.path(out_dir, "deseq_log2FC_by_3prime_category_merged.pdf"), width = 10, height = 4)


#### NEXT: log2FC of first exonic pA- reads, cytoplasmic vs nuclear DHSs

## cytoplasmic: rank below the 4th decile; nuclear: rank above the 7th decile
## (deciles of the first exonic pA- ranks)
tab.lims <- tab.categories %>%
    dplyr::filter(label_3prime_pA %in% c("first_exon_pA-")) %>%
    group_by(nqtile.per.label) %>%
    summarize(m1 = max(overall_rank), m2 = min(overall_rank)) %>%
    dplyr::filter(nqtile.per.label == 4 | nqtile.per.label == 7) %>%
    mutate(m = ifelse(nqtile.per.label == 4, m1, m2)) %>%
    dplyr::select(m) %>%
    mutate(lim.cyto = min(m), lim.nuc = max(m)) %>%
    .[1, c(2, 3)]

deseq.all <- deseq.results$df_lfc_shrunk %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    drop_na() %>%
    mutate(metric = metric %>% str_replace("lfc_", "")) %>%
    separate(DHS, into = c("dhs", "label_3prime_pA"), sep = "\\.") %>%
    dplyr::filter(label_3prime_pA %in% c("gradual_degradation_pA-", "first_exon_pA-")) %>%
    dplyr::select(-lfcSE, -pvalue) %>%
    dplyr::rename(geno = metric, DHS = dhs) %>%
    left_join(tab.categories) %>%
    mutate(label = case_when(
               overall_rank < tab.lims$lim.cyto ~ "cytoplasmic",
               overall_rank > tab.lims$lim.nuc ~ "nuclear",
               TRUE ~ NA_character_
           )) %>%
    mutate(significant = padj < 0.05) %>%
    dplyr::select(-nqtile.per.label) %>%
    dplyr::select(significant, label, geno, DHS, log2FoldChange, label_3prime_pA, baseMean) %>%
    dplyr::filter(geno %in% c("RBM7", "ZCCHC8")) %>%
    drop_na()

g <- ggplot(deseq.all %>% dplyr::filter(grepl("first", label_3prime_pA)), aes(x = label, y = log2FoldChange, fill = geno)) +
    geom_violin(colour = "transparent") +
    theme_bw() +
    coord_cartesian(ylim = c(-2, 2)) +
    ggtitle("Distribution log2FC (NEXT subunit/WT)") +
    ylab("log2 fold change") +
    scale_fill_manual(values = c("#b1cdc3", "#886e97"))
ggsave(file.path(out_dir, "violin_log2FoldChange_NEXT.pdf"), g, width = 5, height = 3)


#### PAXT: DESeq2 MTR4, Z1-FL vs WT, per DHS x category (TLDRSeq_Z1_MTR4_sup)

## chrM scale factors (02); sample names without "-" (mES_Z1-FL_rep1 -> mES_Z1FL_rep1)
df_scale_factor.PAXT <- read.table(file.path(GENIC_DIR, "TLDRSeq_Z1_MTR4_sup", "library_size_factors", "library_size_factors.tsv"), sep = "\t", header = TRUE) %>%
    transmute(scale_factor = scale_factor_chrM,
              sample = sample %>% str_replace("-", "")) %>%
    dplyr::filter(grepl("WT|MTR4|Z1FL", sample))

tab.PAXT <- count_reads_per_category("TLDRSeq_Z1_MTR4_sup") %>%
    dplyr::select(DHS, label_3prime_pA, n, sample) %>%
    dplyr::filter(grepl("WT|MTR4|Z1-FL", sample)) %>%
    group_by(DHS, sample, label_3prime_pA) %>%
    summarize(n = sum(n)) %>%
    mutate(sample = sample %>% str_replace("-", "")) %>%
    mutate(rows = paste0(DHS, ".", label_3prime_pA)) %>%
    reshape2::dcast(rows ~ sample, value.var = "n", fill = 0) %>%
    column_to_rownames("rows") %>%
    dplyr::filter(rowSums(.) > 20)

df.meta.PAXT <- colnames(tab.PAXT) %>%
    data.frame(sample = .) %>%
    mutate(geno = sample %>% str_replace("mES_", "") %>% str_replace("_rep.", "")) %>%
    column_to_rownames("sample") %>%
    mutate(geno = factor(geno, c("WT", "MTR4", "Z1FL")))

deseq.results.PAXT <- run_deseq(tab.PAXT,
                                df.meta.PAXT,
                                c("MTR4", "Z1FL"),
                                df_scale_factor.PAXT)

deseq.results.PAXT.df <- deseq.results.PAXT$df_lfc_raw %>%
    separate(DHS, into = c("garb", "label_3prime_pA"), sep = "\\.") %>%
    mutate(label_3prime = label_3prime_pA %>% str_replace("_pA\\+", "") %>% str_replace("_pA-", "")) %>%
    mutate(pA_status = ifelse(grepl("_pA\\+", label_3prime_pA), "pA", "no_pA")) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    dplyr::select(label_3prime_pA, pA_status, log2FoldChange, geno, garb, baseMean, padj) %>%
    mutate(significant = padj < 0.05) %>%
    dplyr::filter(label_3prime_pA %in% c("first_exon_pA-", "gradual_degradation_pA-")) %>%
    mutate(label_3prime_pA = case_when(
               label_3prime_pA == "first_exon_pA-" ~ "First exonic (pA-)",
               label_3prime_pA == "gradual_degradation_pA-" ~ "Intermediate (pA-)"
           )) %>%
    mutate(label_3prime_pA = factor(label_3prime_pA,
                                    levels = c("Intermediate (pA-)", "First exonic (pA-)")))

g <- ggplot(deseq.results.PAXT.df, aes(x = label_3prime_pA, y = log2FoldChange, fill = geno)) +
    geom_violin(colour = "transparent") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "3' end category", y = "log2 Fold Change (vs WT)") +
    scale_fill_manual(values = c("#E07A5F", "#81B29A"))
ggsave(g, filename = file.path(out_dir, "deseq_log2FC_by_3prime_category_Z1_MTR4.pdf"), width = 5, height = 4)

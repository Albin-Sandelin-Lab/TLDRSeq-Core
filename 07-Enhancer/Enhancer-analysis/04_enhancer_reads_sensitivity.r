## Enhancer DHS reads: sensitivity to RBM7 / ZCCHC8 depletion (DESeq2 per DHS and strand),
## pA fraction, read length and exon number of upregulated enhancers, and log2FC per read
## length category.
##
## Usage (from 07-Enhancer/Enhancer-analysis/): Rscript 04_enhancer_reads_sensitivity.r
##
## Input : $ENHANCER_DIR/coverage_and_chip/DHS_reads_df.tsv (01)
##         SCALE_FACTORS_CSV, CHRM_DEPTH (config/paths.R)
## Output (in $ENHANCER_DIR/reads_sensitivity/):
##   deseq_results_DHS.rds                            DESeq2 RBM7 / ZCCHC8 vs WT per DHS and strand
##   DHS_up_read_exonic_length_density_plot_pA.pdf    median exonic read length of upregulated DHSs
##                                                    (WT: all tested DHSs), pA+ / pA-
##   DHS_up_read_exon_number_fraction_histogram_pA.pdf  median exon number, same DHSs
##   DHS_read_pA_fraction_histogram_same_scale.pdf    pA+ read fraction per DHS (> 5 reads)
##   DHS_list_expressed_reads_over_5.txt              DHSs (and strands) with > 5 reads
##   CI_enhancers.rds                                 k-means log2FC intervals
##   DHS_fraction_pA_vs_logFC_violin_plot_2.pdf       pA+ read fraction per log2FC interval, ANOVA
##   length_category/deseq_results_DHS.rds            DESeq2 per DHS, strand and read length category
##   length_category/DHS_log2FC_violin_plot_length_separated_enhancer.pdf
##                                                    log2FC per read length category
##
## pA+ reads: poly(A) tail > 15 nt. Normalized counts: read count / chrM scale factor, divided by
## the number of replicates. Upregulated: log2FC > 1, padj < 0.05.

library(data.table)
library(tidyverse)
library(classInt)

source("../../config/paths.R")
source("../../lib/tldr_reads.R")
source("../../lib/deseq.R")
source("../../04-Normalization/scale_factors.R")

out_dir <- file.path(ENHANCER_DIR, "reads_sensitivity")
dir.create(file.path(out_dir, "length_category"), showWarnings = FALSE, recursive = TRUE)

reads_file <- file.path(ENHANCER_DIR, "coverage_and_chip", "DHS_reads_df.tsv")


#### chrM scale factors and number of replicates per genotype

df_scale_factor <- read.csv(SCALE_FACTORS_CSV) %>%
    mutate(Rep = sample %>% str_replace("scale_", "")) %>%
    dplyr::select(scale_factor, rep = Rep)

n.lib <- df_scale_factor %>%
    mutate(geno = rep %>% str_replace("mES_", "") %>% str_replace("_rep.", "")) %>%
    mutate(repn = str_replace(rep, ".*_rep", "")) %>%
    dplyr::select(geno, repn) %>%
    group_by(geno) %>%
    slice_max(repn) %>%
    as.data.frame

setDT(df_scale_factor)
setkey(df_scale_factor, "rep")
setDT(n.lib)
setkey(n.lib, "geno")


#### Reads

all_reads_dt <- fread(reads_file)

setnames(all_reads_dt,
         c("V2", "V3", "V5", "V6", "V10", "V11", "V13", "V14"),
         c("start", "end", "pA_len", "strand", "n.exon", "exon.lengths", "DHS", "rep"))

all_reads_dt[, pA_status := ifelse(pA_len > 15, "pA", "no_pA")]
all_reads_dt[, geno := rep %>% str_replace("mES_", "") %>% str_replace("_rep.", "")]
all_reads_dt[, repn_chr := str_replace(rep, ".*_rep", "")]

## normalized read counts per DHS, strand, genotype and pA status (by: extra grouping column)
summarize_counts <- function(reads_dt, by = NULL) {
    keys <- c("rep", "pA_status", "DHS", "strand", "geno", "repn_chr", by)
    counts <- reads_dt[, .N, by = keys]
    setkey(counts, "rep")
    counts <- df_scale_factor[counts]
    counts[, scaled_n := N / scale_factor]
    setkey(counts, "geno")
    counts <- n.lib[counts]
    summarized.pA <- counts[,
        .(
            total = sum(N),
            scaled_total = sum(scaled_n)
        ),
        by = c("pA_status", "geno", "repn", "DHS", "strand", by)
    ]
    ## mean per replicate
    summarized.pA[, repn := as.numeric(repn)]
    summarized.pA[, scaled_total_per_lib := scaled_total / repn]
    summarized.pA[, total_per_lib := total / repn]
    summarized.pA[, repn := NULL]
    summarized.merged <- summarized.pA[,
        .(
            scaled_total_per_lib = sum(scaled_total_per_lib)
        ),
        by = c("DHS", "strand", "geno", by)
    ]
    list(pA = summarized.pA, merged = summarized.merged)
}

summ <- summarize_counts(all_reads_dt)
dhs.reads.number.summarized.pA <- summ$pA
dhs.reads.number.summarized.merged <- summ$merged

## read statistics per DHS, strand, genotype and pA status
all_reads_dt[, pA_length_stat := ifelse(pA_status == "no_pA", 0, pA_len)]
all_reads_dt[, exonic.length := sapply(strsplit(exon.lengths, ",", fixed = TRUE), function(x) sum(as.numeric(x)))]

read_statistics <- function(by) {
    all_reads_dt[,
        .(
            mean.exonic.length = mean(exonic.length),
            median.exonic.length = as.double(median(exonic.length)),
            mean.n.exon = mean(n.exon),
            median.n.exon = as.double(median(n.exon)),
            median.pA.length = as.double(median(pA_length_stat)),
            mean.genomic.length = mean(end - start),
            median.genomic.length = as.double(median(end - start))
        ),
        by = by
    ]
}
dhs.reads.statistics.pA <- read_statistics(c("geno", "pA_status", "DHS", "strand"))
dhs.reads.statistics.merged <- read_statistics(c("geno", "DHS", "strand"))


#### DESeq2 per DHS and strand: RBM7, ZCCHC8 vs WT

dhs.reads.numbers.nopA <- all_reads_dt[, .N, by = .(rep, DHS, strand)]

dhs.prep.FC <- dhs.reads.numbers.nopA %>%
    mutate(DHS.strand = paste(DHS, strand, sep = ".")) %>%
    dplyr::select(DHS.strand, n = N, rep) %>%
    dplyr::filter(!grepl("ZC3H4|ZC3H18|ARS2", rep)) %>%
    reshape2::dcast(DHS.strand ~ rep, value.var = "n", fill = 0) %>%
    column_to_rownames("DHS.strand")

df.meta <- colnames(dhs.prep.FC) %>%
    data.frame(sample = .) %>%
    mutate(geno = sample %>% str_replace("mES_", "") %>% str_replace("_rep.", "")) %>%
    column_to_rownames("sample") %>%
    mutate(geno = factor(geno, c("WT", "RBM7", "ZCCHC8")))

deseq.results <- run_deseq(dhs.prep.FC, df.meta, candidate_geno = c("RBM7", "ZCCHC8"))
saveRDS(deseq.results, file.path(out_dir, "deseq_results_DHS.rds"))

list.dhs.tested <- deseq.results$df_lfc_shrunk$DHS %>% unique

deseq.results.df <- deseq.results$df_lfc_shrunk %>%
    dplyr::select(DHS, log2FoldChange, metric, padj) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    mutate(geno = factor(geno))

## upregulated DHSs (WT: all tested DHSs)
deseq.results.df_up <- deseq.results.df %>%
    drop_na %>%
    dplyr::filter(log2FoldChange > 1 & padj < 0.05) %>%
    dplyr::select(DHS, geno) %>%
    rbind.data.frame(data.frame(DHS = list.dhs.tested, geno = "WT")) %>%
    mutate(status = "upregulated") %>%
    remove_rownames


#### Exonic read length of upregulated DHSs

dhs.up.read.length.pA.df <- dhs.reads.statistics.pA %>%
    mutate(DHS = paste(DHS, strand, sep = ".")) %>%
    dplyr::filter(geno != "ZC3H4") %>%
    mutate(geno = factor(geno, levels = c("WT", "RBM7", "ZCCHC8"))) %>%
    dplyr::select(DHS, geno, median.genomic.length, median.exonic.length, pA_status) %>%
    left_join(deseq.results.df_up, by = c("DHS", "geno")) %>%
    dplyr::filter(!is.na(status))

g <- ggplot(dhs.up.read.length.pA.df, aes(x = median.exonic.length, colour = geno)) +
    geom_density() +
    theme_bw() +
    xlim(0, 3000) +
    facet_grid(rows = vars(pA_status)) +
    scale_colour_manual(values = c('#828282', "#3D405B", "#F2CC8F")) +
    xlab("Median exonic read length") +
    ggtitle("DHS read exonic length distribution")
ggsave(g, filename = file.path(out_dir, "DHS_up_read_exonic_length_density_plot_pA.pdf"), width = 6, height = 4)


#### Exon number of upregulated DHSs

dhs.up.read.exon.pA.df <- dhs.reads.statistics.pA %>%
    dplyr::filter(geno != "ZC3H4") %>%
    mutate(geno = factor(geno, levels = c("WT", "RBM7", "ZCCHC8"))) %>%
    left_join(deseq.results.df_up %>% separate(DHS, into = c("DHS", "strand"), sep = "\\."),
              by = c("DHS", "geno", "strand")) %>%
    dplyr::filter(!is.na(status)) %>%
    dplyr::select(DHS, geno, strand, median.n.exon, pA_status) %>%
    mutate(median.n.exon = round(median.n.exon + 0.1)) %>%
    mutate(group_exon = case_when(
               median.n.exon == 1 ~ "1",
               median.n.exon == 2 ~ "2",
               median.n.exon >= 3 & median.n.exon <= 5 ~ "3-5",
               median.n.exon >= 6 & median.n.exon <= 10 ~ "6-10",
               median.n.exon > 10 ~ ">10"
           )) %>%
    mutate(group_exon = factor(group_exon, levels = c("1", "2", "3-5", "6-10", ">10"))) %>%
    group_by(geno, pA_status) %>%
    add_count() %>%
    ungroup() %>%
    group_by(geno, pA_status, group_exon, n) %>%
    dplyr::summarize(count = n()) %>%
    mutate(fraction = count / n)

g <- ggplot(dhs.up.read.exon.pA.df, aes(x = group_exon, fill = geno, y = fraction)) +
    geom_histogram(stat = "identity", position = position_dodge()) +
    theme_bw() +
    facet_grid(rows = vars(pA_status)) +
    scale_fill_manual(values = c('#828282', "#3D405B", "#F2CC8F")) +
    xlab("Median exon number") +
    ggtitle("DHS read exon number fraction")
ggsave(g, filename = file.path(out_dir, "DHS_up_read_exon_number_fraction_histogram_pA.pdf"), width = 6, height = 4)


#### pA+ read fraction per DHS (> 5 normalized reads)

dhs.pA.fraction.df <- dhs.reads.number.summarized.pA %>%
    dplyr::filter(geno %in% c("WT", "RBM7", "ZCCHC8")) %>%
    dplyr::select(DHS, geno, scaled_total_per_lib, pA_status, strand) %>%
    reshape2::dcast(DHS + geno + strand ~ pA_status, value.var = "scaled_total_per_lib", fill = 0) %>%
    mutate(fraction = pA / (pA + no_pA)) %>%
    left_join(dhs.reads.number.summarized.merged, by = c("DHS", "geno", "strand")) %>%
    dplyr::filter(scaled_total_per_lib > 5) %>%
    mutate(geno = factor(geno, levels = c("WT", "RBM7", "ZCCHC8")))

g <- ggplot(dhs.pA.fraction.df, aes(x = fraction, colour = geno, fill = geno,
                                    y = after_stat(count))) +
    geom_histogram(binwidth = 0.05) +
    theme_bw() +
    facet_grid(row = vars(geno)) +
    scale_colour_manual(values = c('#828282', "#3D405B", "#F2CC8F")) +
    scale_fill_manual(values = c('#828282', "#3D405B", "#F2CC8F")) +
    xlab("pA read fraction") +
    ggtitle("DHS read pA fraction distribution")
ggsave(g, filename = file.path(out_dir, "DHS_read_pA_fraction_histogram_same_scale.pdf"), width = 6, height = 4)


#### pA+ read fraction per log2FC interval (k-means)

dhs.read.length.merged.df <- dhs.reads.statistics.merged %>%
    dplyr::filter(geno %in% c("WT", "RBM7", "ZCCHC8")) %>%
    mutate(geno = factor(geno, levels = c("WT", "RBM7", "ZCCHC8"))) %>%
    dplyr::select(DHS, geno, strand, median.exonic.length)

dhs.list.expressed <- dhs.reads.number.summarized.merged %>%
    dplyr::filter(scaled_total_per_lib > 5) %>%
    dplyr::filter(geno %in% c("WT", "RBM7", "ZCCHC8"))
write.table(dhs.list.expressed %>% dplyr::select(DHS, strand) %>% distinct,
            file = file.path(out_dir, "DHS_list_expressed_reads_over_5.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

logFC.readLength.df <- left_join(dhs.read.length.merged.df,
                                 deseq.results.df %>% separate(DHS, into = c("DHS", "strand"), sep = "\\."),
                                 by = c("DHS", "strand", "geno")) %>%
    left_join(dhs.pA.fraction.df %>% dplyr::select(DHS, strand, geno, fraction),
              by = c("DHS", "strand", "geno")) %>%
    dplyr::select(-padj, -metric) %>%
    mutate(fraction = ifelse(is.na(fraction), 0, fraction)) %>%
    drop_na %>%
    dplyr::rename(fractionpA = fraction) %>%
    left_join(dhs.list.expressed, by = c("DHS", "geno", "strand")) %>%
    drop_na() %>%
    mutate(geno = factor(geno, levels = c("RBM7", "ZCCHC8")))

## 6 log2FC intervals (k-means)
cI <- classIntervals(logFC.readLength.df$log2FoldChange, n = 6, style = "kmeans")
kmeans_breaks <- round(cI$brks, 2)
kmeans_breaks[1] <- kmeans_breaks[1] - 0.01
kmeans_breaks[length(kmeans_breaks)] <- kmeans_breaks[length(kmeans_breaks)] + 0.01
left_bounds <- kmeans_breaks[-length(kmeans_breaks)]
right_bounds <- kmeans_breaks[-1]
interval_strings <- sprintf("[%s, %s)", left_bounds, right_bounds)
num_intervals <- length(interval_strings)
interval_strings[num_intervals] <- sprintf("[%s, %s]", left_bounds[num_intervals], right_bounds[num_intervals])
indices <- findInterval(logFC.readLength.df$log2FoldChange, kmeans_breaks, rightmost.closed = TRUE)
interval_assignments <- factor(interval_strings[indices], levels = interval_strings)
saveRDS(cI, file.path(out_dir, "CI_enhancers.rds"))

logFC.readLength.df <- logFC.readLength.df %>%
    mutate(id3 = interval_assignments)

## ANOVA: pA+ fraction ~ log2FC interval, per genotype
ano.list <- list()
for (geno.elt in unique(logFC.readLength.df$geno))
{
    ano.list.df <- logFC.readLength.df %>%
        dplyr::filter(geno.elt == geno)
    ano.list[[geno.elt]] <- lm(fractionpA ~ id3, data = ano.list.df) %>% anova
}

ano.df <- ano.list %>% sapply(getElement, "Pr(>F)") %>% .[1, ] %>% data.frame(geno = names(.), pval = .) %>%
    mutate(cat = paste("ANOVA pval =", formatC(pval, format = "e", digits = 2))) %>%
    remove_rownames() %>%
    dplyr::select(-pval)

p1 <- ggplot(logFC.readLength.df, aes(y = fractionpA, x = id3, fill = geno)) +
    geom_violin(color = "black", size = 0.2) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.placement = "outside") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    stat_summary(fun = mean, geom = "point", shape = 1, size = 2, color = "black") +
    stat_summary(fun = median, geom = "point", shape = 2, size = 2, color = "black") +
    facet_grid(row = vars(geno), switch = "x") +
    scale_fill_manual(values = c("#3D405B", "#F2CC8F")) +
    ylab("fraction pA reads") +
    ylim(c(0, 1)) +
    xlab("log2FC (in k-mean intervals)") +
    ## ANOVA p-value written above the second to last interval
    geom_text(aes(label = cat), data = logFC.readLength.df %>% left_join(ano.df) %>% distinct(cat, .keep_all = TRUE),
              size = 3, y = 0.95, x = interval_strings[num_intervals - 1])
ggsave(p1, filename = file.path(out_dir, "DHS_fraction_pA_vs_logFC_violin_plot_2.pdf"), width = 6, height = 4)


#### log2FC per read length category

all_reads_dt[, group_cat := fcase(
                   end - start < 100, "0-100",
                   end - start < 150, "100-150",
                   end - start < 200, "150-200",
                   end - start < 250, "200-250",
                   end - start < 300, "250-300",
                   end - start < 350, "300-350",
                   end - start < 400, "350-400",
                   end - start < 500, "400-500",
                   end - start < 600, "500-600",
                   end - start < 700, "600-700",
                   end - start < 800, "700-800",
                   end - start < 1000, "800-1000",
                   end - start < 1500, "1000-1500",
                   default = ">1500")
             ]

group_cat.levels <- c("0-100", "100-150", "150-200", "200-250", "250-300", "300-350", "350-400",
                      "400-500", "500-600", "600-700", "700-800", "800-1000", "1000-1500", ">1500")

summ.cat <- summarize_counts(all_reads_dt, by = "group_cat")

## DESeq2 per DHS, read length category and strand
dhs.prep.FC.cat <- all_reads_dt[, .N, by = .(rep, DHS, strand, group_cat)] %>%
    dplyr::select(DHS, strand, group_cat, n = N, rep) %>%
    mutate(DHS.group_cat = paste0(DHS, ".", group_cat, ".", strand)) %>%
    dplyr::filter(!grepl("ZC3H4|ZC3H18|ARS2", rep)) %>%
    reshape2::dcast(DHS.group_cat ~ rep, value.var = "n", fill = 0) %>%
    column_to_rownames("DHS.group_cat")

df.meta.cat <- colnames(dhs.prep.FC.cat) %>%
    data.frame(sample = .) %>%
    mutate(geno = sample %>% str_replace("mES_", "") %>% str_replace("_rep.", "")) %>%
    column_to_rownames("sample") %>%
    mutate(geno = factor(geno, c("WT", "RBM7", "ZCCHC8")))

deseq.results.cat <- run_deseq(dhs.prep.FC.cat, df.meta.cat, candidate_geno = c("RBM7", "ZCCHC8"))
saveRDS(deseq.results.cat, file.path(out_dir, "length_category", "deseq_results_DHS.rds"))

deseq.results.cat.df <- deseq.results.cat$df_lfc_shrunk %>%
    dplyr::select(DHS, log2FoldChange, metric, padj) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    mutate(geno = factor(geno, levels = c("RBM7", "ZCCHC8")))

## DHS x category x strand with > 2 normalized pA+ or pA- reads
deseq.results.df.joined <- deseq.results.cat.df %>%
    separate(DHS, into = c("DHS", "group_cat", "strand"), sep = "\\.") %>%
    left_join(summ.cat$merged,
              by = c("DHS", "geno", "strand", "group_cat")) %>%
    dplyr::filter(!is.na(scaled_total_per_lib)) %>%
    left_join(summ.cat$pA %>% dplyr::select(-total_per_lib, -total, -scaled_total),
              by = c("DHS", "geno", "group_cat", "strand"),
              multiple = "all") %>%
    dplyr::rename(tpm.merged = scaled_total_per_lib.x,
                  tpm.pA = scaled_total_per_lib.y) %>%
    dplyr::filter(tpm.pA > 2) %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate(status = ifelse(padj < 0.05 & log2FoldChange > 1, "upregulated", "unchanged"))

deseq.results.df.joined.plot.violin <- deseq.results.df.joined %>%
    dplyr::select(-status, -tpm.pA, -tpm.merged, -pA_status) %>%
    distinct() %>%
    mutate(group_cat = factor(group_cat, levels = group_cat.levels)) %>%
    mutate(geno = factor(geno, levels = c("RBM7", "ZCCHC8")))

g <- ggplot(deseq.results.df.joined.plot.violin, aes(x = group_cat, y = log2FoldChange, fill = geno)) +
    geom_violin(color = "black", size = 0.2) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.placement = "outside") +
    stat_summary(fun = mean, geom = "point", shape = 1, size = 2, color = "black") +
    stat_summary(fun = median, geom = "point", shape = 2, size = 2, color = "black") +
    facet_grid(row = vars(geno), switch = "x") +
    scale_fill_manual(values = c("#3D405B", "#F2CC8F")) +
    xlab("DHS length category") +
    ylab("log2 Fold Change vs WT") +
    ggtitle(paste0("DHS log2FC upon exosome related factor KD - enhancer reads"))
ggsave(g, filename = file.path(out_dir, "length_category", "DHS_log2FC_violin_plot_length_separated_enhancer.pdf"), width = 10, height = 4)

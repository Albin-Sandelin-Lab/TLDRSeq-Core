## Library size (mapped reads) of each TLDRSeq_sup replicate (WT, RBM7, ZCCHC8).
##
## Usage (from 08-QC/): Rscript 01_library_size.r
##
## Input : $RESULTS_DIR/TLDRSeq_sup/concatenate_per_replicate/<sample>.bam(.bai) (01-TLDRSeq-alignment)
## Output (in $QC_DIR/library_size/):
##   lib_size.tsv               BAM file, number of mapped reads (from the BAM index, = samtools idxstats)
##   library_size_barplot.pdf   library size per replicate

library(tidyverse)
library(Rsamtools)

source("../config/paths.R")

in_dir  <- file.path(RESULTS_DIR, "TLDRSeq_sup", "concatenate_per_replicate")
out_dir <- file.path(QC_DIR, "library_size")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#### Mapped reads per BAM

bam.files <- list.files(in_dir, pattern = "\\.bam$")
lib.size <- sapply(bam.files, function(f) sum(idxstatsBam(file.path(in_dir, f))$mapped))

write.table(data.frame(lib = bam.files, n = lib.size), file.path(out_dir, "lib_size.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#### Plot

lib.size.df <- read.csv(file.path(out_dir, "lib_size.tsv"), header = FALSE, sep = "\t") %>%
    setNames(c("lib", "n")) %>%
    dplyr::filter(grepl("RBM7|ZCCHC8|WT", lib)) %>%
    separate(lib, into = c("mes", "geno", "rep", "bam"), sep = "[\\._]") %>%
    dplyr::select(-mes, -bam) %>%
    mutate(rep = str_replace(rep, "rep", "")) %>%
    arrange(geno, rep) %>%
    mutate(geno = factor(geno, levels = c("WT", "RBM7", "ZCCHC8")))

g <- ggplot(lib.size.df, aes(x = geno, y = n, colour = geno)) +
    geom_point() +
    scale_colour_manual(values = c('#828282', "#3D405B", "#F2CC8F")) +
    scale_x_discrete(drop = TRUE) +
    labs(x = "Replicate", y = "Library Size (reads)", fill = "Genotype") +
    theme_bw() +
    ylim(5.5e7, 8e7) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(g, filename = file.path(out_dir, "library_size_barplot.pdf"), width = 4, height = 4)

## Sequence logos around the TSSs of the enhancer transcripts, per TSS distance quantile.
##
## Usage (from 07-Enhancer/Enhancer-analysis/): Rscript 05_enhancer_tss_logo.r
##
## Sequences from -100 to +20 nt around the + strand TSS and the - strand TSS of each enhancer
## DHS (on the transcript strand; sequences with N excluded), pooled per TSS distance quantile.
## Logo height: nucleotide frequency (pseudocount 0.01) x information content.
##
## Input : $ENHANCER_DIR/coverage_and_chip/{tss_distance.tsv, DHS_ntile_distance.tsv} (01)
##         GENOME_FA (config/paths.R)
## Output (in $ENHANCER_DIR/tss_logo/):
##   logo2_all_ntiles_all.pdf  one logo per TSS distance quantile

library(GenomicRanges)
library(Biostrings)
library(parallel)
library(tidyverse)

source("../../config/paths.R")

out_dir <- file.path(ENHANCER_DIR, "tss_logo")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genome <- Biostrings::readDNAStringSet(GENOME_FA, format = "fasta", use.names = TRUE)

dhs.order.for.center <- read.table(file.path(ENHANCER_DIR, "coverage_and_chip", "DHS_ntile_distance.tsv"), header = TRUE, sep = "\t")
dhs_center.df <- read.table(file.path(ENHANCER_DIR, "coverage_and_chip", "tss_distance.tsv"), sep = "\t", header = TRUE)


#### Sequences around the + and - strand TSSs (+/- 100 nt)

TSS.plus <- dhs_center.df %>%
    dplyr::select(DHS, TSS_plus) %>%
    left_join(dhs.order.for.center %>% dplyr::select(DHS, ntile.dist), by = "DHS") %>%
    mutate(seqnames = DHS %>% as.character() %>% strsplit("_", fixed = TRUE) %>% sapply("[", 1),
           start = TSS_plus - 100,
           end = TSS_plus + 100,
           strand = "+") %>%
    dplyr::select(-DHS, -TSS_plus) %>%
    plyranges::as_granges()

TSS.minus <- dhs_center.df %>%
    dplyr::select(DHS, TSS_minus) %>%
    left_join(dhs.order.for.center %>% dplyr::select(DHS, ntile.dist), by = "DHS") %>%
    mutate(seqnames = DHS %>% as.character() %>% strsplit("_", fixed = TRUE) %>% sapply("[", 1),
           start = TSS_minus - 100,
           end = TSS_minus + 100,
           strand = "-") %>%
    dplyr::select(-DHS, -TSS_minus) %>%
    plyranges::as_granges()

TSS.all <- c(TSS.plus, TSS.minus) %>%
    split(f = .$ntile.dist)

sequences.list <- lapply(TSS.all, function(x, y) BSgenome::getSeq(y, x), genome) %>%
    lapply(function(x) x[!grepl("[^ACTG]", as.character(x))])


#### Nucleotide counts per position

## counts of A, C, G, T per position; positions trimmed by reduce.left / reduce.right
nt.freq <- function(seq.logo, reduce.left = 0, reduce.right = 0)
{
    freq <- seq.logo %>%
        as.character() %>%
        str_split(., "") %>%
        do.call(rbind, .) %>%
        as.data.frame %>%
        as.list() %>%
        lapply(paste, collapse = "") %>%
        lapply(function(x) sapply(list(A = "A", C = "C", G = "G", T = "T"), function(y) str_count(x, y))) %>%
        do.call(cbind.data.frame, .) %>%
        t %>%
        as.data.frame %>%
        remove_rownames %>%
        mutate(position = 1:nrow(.)) %>%
        dplyr::filter(position %in% seq(reduce.left, max(position) - reduce.right)) %>%
        mutate(position = position - reduce.left) %>%
        reshape2::melt(id.vars = c("position"))
    return(freq)
}

## -100 to +20 nt
freq.logo <- mclapply(sequences.list, nt.freq, 0, 80, mc.cores = 10)

all.freq <- mapply(mutate, freq.logo,
                   seq_group = names(freq.logo), SIMPLIFY = FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    dplyr::rename(nt = variable) %>%
    group_by(seq_group, position) %>%
    mutate(n.seq = sum(value)) %>%
    group_by(nt, seq_group) %>%
    mutate(nt.abundance = sum(value)) %>%
    ungroup %>%
    group_by(seq_group) %>%
    mutate(nt.tot = sum(value)) %>%
    ungroup %>%
    mutate(freq = value / n.seq, freq.norm = value / (n.seq * (nt.abundance / nt.tot)))


#### Logos

## frequency matrices (pseudocount 0.01), one per quantile
all.logo2.towrite <- all.freq %>%
    dplyr::select(seq_group, position, nt, freq) %>%
    group_by(position, seq_group) %>%
    mutate(freq = freq + 0.01) %>%
    mutate(tot = sum(freq), freq = (freq / tot)) %>%
    ungroup %>%
    reshape2::dcast(seq_group + position ~ nt, value.var = "freq") %>%
    remove_rownames %>%
    split(f = .$seq_group) %>%
    lapply(dplyr::select, A, C, G, T) %>%
    lapply(remove_rownames)

## letter height = frequency x information content (2 - entropy)
p_list2 <- all.logo2.towrite %>%
    lapply(function(x) t(as.matrix(x) * matrix(rep(2 + rowSums((as.matrix(x)) * log(as.matrix(x), 2)), 4), ncol = 4)))

g <- ggseqlogo::ggseqlogo(p_list2 %>% setNames(paste0("quantile ", names(.))),
                          method = "custom", seq_type = "dna") +
    facet_wrap(~seq_group, ncol = 1, dir = "v") +
    scale_x_continuous(breaks = seq(1, 121, 10), labels = seq(-100, 20, 10), guide = guide_axis(n.dodge = 1))
ggpubr::ggexport(g, filename = file.path(out_dir, "logo2_all_ntiles_all.pdf"), width = 10, height = 15)

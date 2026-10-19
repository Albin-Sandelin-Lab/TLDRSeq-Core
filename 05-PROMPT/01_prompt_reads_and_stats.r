## Collect the TLDR-Seq reads of each DHS, keep PROMPT reads and summarize them.
##
## Usage (from 05-PROMPT/): Rscript 01_prompt_reads_and_stats.r
##
## DHSs: uniquely associated with one RefSeq select gene TSS, not at a TTS, not overlapping
## an MGI ncRNA, and with a non-ambiguous gene strand.
## PROMPT reads: reads of these DHSs on the strand opposite to the gene.
##
## Output (in $PROMPT_DIR/prompt_reads/):
##   DHS_all_reads_raw.tsv                   all reads of the selected DHSs (both strands), used by 03 and 05
##   DHS_reads_counts_nopA.tsv               PROMPT read counts per DHS, strand and replicate (DESeq2 input)
##   DHS_reads_number_summarized_pA.tsv      normalized PROMPT read counts per DHS, genotype and pA status
##   DHS_reads_number_summarized_merged.tsv  same, pA+ and pA- summed
##   DHS_reads_statistics_pA.tsv             read length / exon number / pA length per DHS, genotype, pA status
##   DHS_reads_statistics_merged.tsv         same, pA+ and pA- together
##
## pA+ reads: poly(A) tail > 15 nt.
## Normalized counts: read count / chrM scale factor, divided by the number of replicates.

library(rtracklayer)
library(GenomicRanges)
library(parallel)
library(data.table)
library(tidyverse)

source("../config/paths.R")
source("../lib/tldr_reads.R")

out_dir <- file.path(PROMPT_DIR, "prompt_reads")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


#### Select DHSs

df_DHS_prompt_map <- readRDS(DHS_ANNOTATION_RDS) %>%
    filter(filter_refseq_select_gene_TSS == "UNIQ", # DHS uniquely mapped to a gene TSS
           filter_refseq_select_gene_TTS == "NONE" |
           (filter_refseq_select_gene_TTS == "UNIQ" & annots_refseq_select_gene_TTS == annots_refseq_select_gene_TSS), # DHS not mapped to a TTS
           filter_refseq_select_gene_full %in% c("UNIQ", "ANNOT-MULTI_DHS"),
           filter_gencode_gene_full %in% c("UNIQ", "ANNOT-MULTI_DHS"),
           filter_gencode_gene_full_vs_DHS2000 %in% c("UNIQ", "ANNOT-MULTI_DHS"),
           filter_mgi_ncRNA == "NONE" # DHS not mapped to a ncRNA
           )

dhs.list <- df_DHS_prompt_map$DHS %>% unique

dhs.gene.strand <- read.table(DHS_GENE_STRAND, sep = "\t", header = TRUE)


#### Scale factors (chrM) and number of replicates per genotype

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


#### All reads of the selected DHSs

get.dhs.reads <- function(dhs, geno = ".") {
    tmp2 <- tmp_replicates[grepl(geno, tmp_replicates)]
    DHS_reads <- read_TLDR_reads_by_list_of_DHS(query_list_DHS = dhs, query_list_replicates = tmp2)
    if (nrow(DHS_reads) == 0) {
        return(NULL)
    }
    return(DHS_reads)
}

dhs.read.list <- mclapply(dhs.list, get.dhs.reads, mc.cores = 40)
dhs.read.list <- dhs.read.list[!sapply(dhs.read.list, is.null)]

all_reads_dt <- rbindlist(dhs.read.list, use.names = FALSE)
fwrite(all_reads_dt, file.path(out_dir, "DHS_all_reads_raw.tsv"), sep = "\t")
rm(dhs.read.list)
gc()

setnames(all_reads_dt,
         c("V2", "V3", "V5", "V6", "V10", "V11", "V13", "V14"),
         c("start", "end", "pA_len", "strand", "n.exon", "exon.lengths", "DHS", "rep"))


#### PROMPT reads: strand opposite to the gene

dhs.select <- dhs.gene.strand %>%
    dplyr::filter(gene_strand != "ambiguous") %>%
    .$DHS %>%
    unique

all_reads_dt[, dhs.strand := paste0(DHS, ".", strand)]
all_reads_dt <- all_reads_dt[DHS %in% dhs.select, ]

## DHS.strand = DHS + gene strand: reads on this strand are gene reads
DHS.strand <- paste(dhs.gene.strand$DHS, dhs.gene.strand$gene_strand, sep = ".") %>%
    .[!grepl("ambiguous", .)] %>%
    .[!grepl("\\.NA$", .)]

all_reads_dt <- all_reads_dt[!dhs.strand %in% DHS.strand, ]


#### Read counts per DHS, genotype and pA status

setDT(df_scale_factor)
setkey(df_scale_factor, "rep")
setDT(n.lib)
setkey(n.lib, "geno")

all_reads_dt[, pA_status := ifelse(pA_len > 15, "pA", "no_pA")]
all_reads_dt[, geno := rep %>% str_replace("mES_", "") %>% str_replace("_rep.", "")]
all_reads_dt[, repn_chr := str_replace(rep, ".*_rep", "")]

dhs.reads.counts.pA <- all_reads_dt[, .N, by = .(rep, pA_status, DHS, strand, geno, repn_chr)]

## Normalize by the chrM scale factor
setkey(dhs.reads.counts.pA, "rep")
dhs.reads.counts.pA <- df_scale_factor[dhs.reads.counts.pA]
dhs.reads.counts.pA[, scaled_n := N / scale_factor]

## Add the number of replicates of each genotype (repn)
setkey(dhs.reads.counts.pA, "geno")
dhs.reads.counts.pA <- n.lib[dhs.reads.counts.pA]

dhs.reads.number.summarized.pA <- dhs.reads.counts.pA[,
    .(
        total = sum(N),
        scaled_total = sum(scaled_n)
    ),
    by = .(pA_status, geno, repn, DHS, strand)
]

## Mean per replicate
dhs.reads.number.summarized.pA[, repn := as.numeric(repn)]
dhs.reads.number.summarized.pA[, scaled_total_per_lib := scaled_total / repn]
dhs.reads.number.summarized.pA[, total_per_lib := total / repn]
dhs.reads.number.summarized.pA[, repn := NULL]

dhs.reads.number.summarized.merged <- dhs.reads.number.summarized.pA[,
    .(
        scaled_total_per_lib = sum(scaled_total_per_lib)
    ),
    by = .(DHS, strand, geno)
]


#### Read statistics per DHS, genotype and pA status

all_reads_dt[, pA_length_stat := ifelse(pA_status == "no_pA", 0, pA_len)]
all_reads_dt[, exonic.length := sapply(strsplit(exon.lengths, ",", fixed = TRUE), function(x) sum(as.numeric(x)))]

dhs.reads.statistics.pA <- all_reads_dt[,
    .(
        mean.exonic.length = mean(exonic.length),
        median.exonic.length = as.double(median(exonic.length)),
        mean.n.exon = mean(n.exon),
        median.n.exon = as.double(median(n.exon)),
        median.pA.length = as.double(median(pA_length_stat)),
        mean.genomic.length = mean(end - start),
        median.genomic.length = as.double(median(end - start))
    ),
    by = .(geno, pA_status, DHS, strand)
]

dhs.reads.statistics.merged <- all_reads_dt[,
    .(
        mean.exonic.length = mean(exonic.length),
        median.exonic.length = as.double(median(exonic.length)),
        mean.n.exon = mean(n.exon),
        median.n.exon = as.double(median(n.exon)),
        median.pA.length = as.double(median(pA_length_stat)),
        mean.genomic.length = mean(end - start),
        median.genomic.length = as.double(median(end - start))
    ),
    by = .(geno, DHS, strand)
]


#### Raw read counts per replicate (DESeq2 input)

dhs.reads.counts.nopA <- all_reads_dt[, .N, by = .(rep, DHS, strand)]


#### Write

fwrite(dhs.reads.counts.nopA, file.path(out_dir, "DHS_reads_counts_nopA.tsv"), sep = "\t")
fwrite(dhs.reads.number.summarized.pA, file.path(out_dir, "DHS_reads_number_summarized_pA.tsv"), sep = "\t")
fwrite(dhs.reads.statistics.pA, file.path(out_dir, "DHS_reads_statistics_pA.tsv"), sep = "\t")
fwrite(dhs.reads.number.summarized.merged, file.path(out_dir, "DHS_reads_number_summarized_merged.tsv"), sep = "\t")
fwrite(dhs.reads.statistics.merged, file.path(out_dir, "DHS_reads_statistics_merged.tsv"), sep = "\t")

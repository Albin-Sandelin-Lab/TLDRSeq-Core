## PROMPT vs PMT vs first exon vs gene reads, by read length category.
##
## Usage (from 05-PROMPT/): Rscript 03_prompt_vs_pmt_length.r
##
## Read types (per DHS, gene from RefSeq select):
##   prompt : reads on the strand opposite to the gene
##   exon1  : gene-strand reads whose 3' end is in the first exon (minus its last 10 nt)
##   pmt    : gene-strand reads whose 3' end is in the first intron
##   gene   : gene-strand reads whose 3' end is in the gene (+/- 500 bp)
## Length categories: genomic read length (end - start), from 0-100 to >1500.
##
## Input : $PROMPT_DIR/prompt_reads/DHS_all_reads_raw.tsv (01_prompt_reads_and_stats.r)
## Output: $PROMPT_DIR/prompt_vs_pmt_length/
##   DHS_reads_counts_nopA.tsv, DHS_reads_number_summarized_{pA,merged}.tsv,
##   DHS_reads_statistics_{pA,merged}.tsv, deseq_results_DHS.rds, figures
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(parallel)
library(pheatmap)
library(data.table)
library(org.Mm.eg.db)
library(tidyverse)
library(RColorBrewer)
library(GenomicFeatures)
library(IRanges)

source("../config/paths.R")
source("../lib/tldr_reads.R")
source("../lib/deseq.R")
source("../04-Normalization/scale_factors.R")

out_dir <- file.path(PROMPT_DIR, "prompt_vs_pmt_length")
dir.create(out_dir, showWarnings=FALSE, recursive =TRUE)

file_df_DHS_prompt_mapped <- DHS_ANNOTATION_RDS

df_DHS_prompt_map <- readRDS(file_df_DHS_prompt_mapped) %>%
    filter(filter_refseq_select_gene_TSS == "UNIQ", # DHS uniquely mapped to gencode
         filter_refseq_select_gene_TTS == "NONE" |
         (filter_refseq_select_gene_TTS == "UNIQ" & annots_refseq_select_gene_TTS == annots_refseq_select_gene_TSS), # DHS not map to TTS
         filter_refseq_select_gene_full %in% c("UNIQ", "ANNOT-MULTI_DHS"),
         filter_gencode_gene_full %in% c("UNIQ", "ANNOT-MULTI_DHS"),
         filter_gencode_gene_full_vs_DHS2000 %in% c("UNIQ", "ANNOT-MULTI_DHS"),
         filter_mgi_ncRNA == "NONE"# DHS not map to ncRNA
         )

dhs.list <- df_DHS_prompt_map$DHS %>% unique

dhs.gene.strand <- read.table(DHS_GENE_STRAND, sep="\t", header=TRUE)

###################" running some tests
txdb = makeTxDbFromGFF(REFSEQ_SELECT_GTF)

raw_exons_df <- exons(txdb, columns = c("GENEID", "exon_rank")) %>%
    as.data.frame()

clean_exons_df <- raw_exons_df %>%
  tidyr::unnest(cols = c(GENEID, exon_rank))

gr_tss <- clean_exons_df %>%
    filter(exon_rank == 1) %>%
    plyranges::as_granges() %>%
    resize(width = 1, fix = "start")

gr_exon1 <- clean_exons_df %>%
    group_by(GENEID) %>%
    filter(n() >= 2)  %>%
    ungroup()  %>% 
    filter(exon_rank == 1) %>%
    plyranges::as_granges() %>%
    plyranges::filter(width > 50) %>% 
    plyranges::anchor_5p() %>%
    mutate(width = width - 10) %>%
    as.data.frame %>%
    plyranges::as_granges()

### 2. Filter for First intron
gr_intron1 <- clean_exons_df %>%
  filter(exon_rank %in% c(1, 2)) %>%
    group_by(GENEID) %>%
    filter(n() == 2) %>%
    summarise(
    seqnames = unique(seqnames),
    strand   = unique(strand),
    intron_start    = min(end) + 1,    # End of the earlier exon + 1
    intron_end      = max(start) - 1   # Start of the later exon - 1
    ) %>%
    dplyr::select(seqnames, start = intron_start, end = intron_end, strand, GENEID) %>%
    plyranges::as_granges() 

gr_genes <- clean_exons_df %>%
    group_by(GENEID) %>%
    summarise(
        seqnames = unique(seqnames),
        strand   = unique(strand),
        gene_start    = min(start) -500,  
        gene_end      = max(end) + 500 ,    
    ) %>%
    dplyr::select(seqnames, start = gene_start, end = gene_end, strand, GENEID) %>%
    plyranges::as_granges() 

    

######## load scale factor from mitochondria
df_scale_factor = read.csv(SCALE_FACTORS_CSV) %>%
    mutate(Rep = sample %>% str_replace("scale_","")) %>%
    dplyr::select(scale_factor, rep = Rep)

n.lib <- df_scale_factor %>%
    mutate(geno = rep %>% str_replace("mES_","") %>% str_replace("_rep.",""))  %>%
    mutate(repn = str_replace(rep, ".*_rep",""))  %>%
    dplyr::select(geno, repn) %>%
    group_by(geno) %>%
    slice_max(repn) %>%
    as.data.frame

all_reads_dt <- fread(file.path(PROMPT_DIR, "prompt_reads", "DHS_all_reads_raw.tsv"))

setnames(all_reads_dt, 
         c("V2", "V3", "V5", "V6", "V10", "V11", "V13", "V14"), 
         c("start", "end", "pA_len", "strand", "n.exon", "exon.lengths", "DHS", "rep"))

dhs.select <- dhs.gene.strand %>%
    dplyr::filter(gene_strand != "ambiguous") %>%
    .$DHS %>%
    unique

all_reads_dt2 <- all_reads_dt[, dhs.strand := paste0(DHS, ".",strand)]
all_reads_dt2 <- all_reads_dt2[DHS %in% dhs.select,]
DHS.strand <- paste(dhs.gene.strand$DHS, dhs.gene.strand$gene_strand, sep=".") %>%
    .[!grepl("ambiguous", .)] %>% 
    .[!grepl("\\.NA$", .)] 

all_prompt_reads_dt <- all_reads_dt2[! dhs.strand %in% DHS.strand,]
all_gene_reads_dt <- all_reads_dt2[ dhs.strand %in% DHS.strand,]

rm(all_reads_dt)
gc()

### Here, we are going to select the PMT reads only.

all_gene_reads.dhs.gr <- all_gene_reads_dt %>%
    dplyr::select(dhs.strand) %>%
    separate(dhs.strand, into = c("seqnames", "start", "end", "strand"), sep = "[\\._-]") %>%
    mutate(strand = ifelse(strand == "", "-", strand))  %>%
    mutate(start = as.integer(start), end = as.integer(end)) %>%
    plyranges::as_granges()

reads_gr.end <- GRanges(
    seqnames = all_gene_reads_dt$V1, ## Assuming V1 is chromosome
    ranges = IRanges(start = all_gene_reads_dt$start, 
                     end = all_gene_reads_dt$end),
    strand = all_gene_reads_dt$strand
    ) %>%
    plyranges::anchor_3p() %>%
    mutate(width = 1)

##### dhs.df
dhs.gr <- DHS.strand  %>%
    data.frame(DHS = .) %>%
    separate( DHS, into = c("seqnames", "start", "end", "strand"), sep = "[\\._-]") %>%
    mutate(strand = ifelse(strand == "", "-", strand)) %>%
    mutate(start = as.integer(start), end = as.integer(end)) %>%
    plyranges::as_granges() 

refseq.dhs.ov <- findOverlaps(dhs.gr, gr_tss, ignore.strand = FALSE)

dhs.gene.df <- data.frame(
    DHS = DHS.strand[queryHits(refseq.dhs.ov)] %>% str_replace("\\..*$",""),
    refseq_gene = gr_tss$GENEID[subjectHits(refseq.dhs.ov)]
) %>%
    distinct(DHS, .keep_all = TRUE)

#####

ov.read.end.1st.exon <- findOverlaps(reads_gr.end, gr_exon1, ignore.strand = FALSE)
ov.read.end.1st.intron <- findOverlaps(reads_gr.end, gr_intron1, ignore.strand = FALSE)
ov.read.end.genes <- findOverlaps(reads_gr.end, gr_genes, ignore.strand = FALSE)

index.reads.intron <- ov.read.end.1st.intron %>% queryHits
index.reads.exon <- ov.read.end.1st.exon %>% queryHits
index.reads.genes <- ov.read.end.genes %>% queryHits

all_pmt_reads_dt <- all_gene_reads_dt[index.reads.intron, ]
all_exon1_reads_dt <- all_gene_reads_dt[index.reads.exon, ]
all_genes_reads_dt2 <- all_gene_reads_dt[index.reads.genes, ]

name.map.intron <- gr_intron1[subjectHits(ov.read.end.1st.intron)]$GENEID
name.map.exon1 <- gr_exon1[subjectHits(ov.read.end.1st.exon)]$GENEID
name.map.genes <- gr_genes[subjectHits(ov.read.end.genes)]$GENEID

all_pmt_reads_dt[, tts.name := name.map.intron ]
all_exon1_reads_dt[, tts.name := name.map.exon1 ]
all_genes_reads_dt2[, tts.name := name.map.genes ]

all_pmt_reads_dt <-  merge(all_pmt_reads_dt, dhs.gene.df, by="DHS")
all_exon1_reads_dt <-  merge(all_exon1_reads_dt, dhs.gene.df, by="DHS")
all_genes_reads_dt2 <-  merge(all_genes_reads_dt2, dhs.gene.df, by="DHS")

all_pmt_reads_dt <-  all_pmt_reads_dt[refseq_gene ==  tts.name,]
all_exon1_reads_dt <-  all_exon1_reads_dt[refseq_gene ==  tts.name,]
all_genes_reads_dt2 <-  all_genes_reads_dt2[refseq_gene ==  tts.name,]

all_pmt_reads_dt[, dhs.type := paste0(DHS, ".","pmt")]
all_exon1_reads_dt[, dhs.type := paste0(DHS, ".","exon1")]
all_genes_reads_dt2[, dhs.type := paste0(DHS, ".","gene")]
all_prompt_reads_dt[, `:=` (dhs.type = paste0(DHS, ".","prompt"), tts.name = NA, refseq_gene = NA)]

all_reads_dt <- rbind(all_prompt_reads_dt, all_pmt_reads_dt, all_genes_reads_dt2, all_exon1_reads_dt)

all_reads_dt[, group_cat:= fcase(
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

print("Master data.table created.")

setDT(df_scale_factor)
setkey(df_scale_factor, "rep")
setDT(n.lib)
setkey(n.lib, "geno")
###### 3. Read counts per DHS, strand and genotype (pA+ / pA- separated)
print("Calculating read counts and summaries...")

## Add helper columns for pA status and genotype
all_reads_dt[, pA_status := ifelse(pA_len > 15, "pA", "no_pA")]

## Use the same logic  used for n.lib to create a matching 'geno' column
all_reads_dt[, geno := rep %>% str_replace("mES_","") %>% str_replace("_rep.","")]

all_reads_dt[, repn_chr := str_replace(rep, ".*_rep","")]

## read counts per replicate
## It counts reads per group, equivalent to get.read.N.per.rep
dhs.reads.counts.pA <- all_reads_dt[, .N, by = .(rep, pA_status, dhs.type, strand, geno, repn_chr, group_cat)]

setkey(dhs.reads.counts.pA, "rep")
dhs.reads.counts.pA <- df_scale_factor[dhs.reads.counts.pA] ## Join
dhs.reads.counts.pA[, scaled_n := N / scale_factor]

## summarize per genotype (no minimum number of replicates)
## We join with n.lib to get the max rep number
## This join will now work and add the 'repn' column
setkey(dhs.reads.counts.pA, "geno")
dhs.reads.counts.pA <- n.lib[dhs.reads.counts.pA] 

## The rest of the code should work as intended
dhs.reads.counts.pA[, n.samples := .N, by = .(geno, pA_status, strand, dhs.type, group_cat)]
## Filter for DHSs present in all reps
## dhs.reads.counts.filtered <- dhs.reads.counts.pA[n.samples == repn]
dhs.reads.counts.filtered <- dhs.reads.counts.pA
dhs.list.filtered <- dhs.reads.counts.filtered$dhs.type %>% unique
## Now summarize, replacing 'dhs.reads.number.summarized.pA'

dhs.reads.number.summarized.pA <- dhs.reads.counts.filtered[,
    .(
        total = sum(N),
        scaled_total = sum(scaled_n)
    ),
    by = .(pA_status, geno, repn, dhs.type, strand, group_cat) ## Group by all keys
]

                                        ## Calculate per-library totals
dhs.reads.number.summarized.pA[, repn := as.numeric(repn)]
dhs.reads.number.summarized.pA[, scaled_total_per_lib := scaled_total / repn]
dhs.reads.number.summarized.pA[, total_per_lib := total / repn]

dhs.reads.number.summarized.pA[, repn := NULL] ## Drop repn

print("pA read summaries complete.")

dhs.reads.number.summarized.merged <- dhs.reads.number.summarized.pA[,
    .(
        # Sum the 'scaled_total_per_lib' for the pA/no_pA groups
        scaled_total_per_lib = sum(scaled_total_per_lib)
    ),
    by = .(dhs.type, geno, strand, group_cat) # Group by DHS, genotype, and strand
]

###### 4. Read statistics (length, exon number, pA length)
print("Calculating read statistics...")

# Calculate exonic.length.
# The str_split/sapply is unavoidable, but doing it once on the
# master table is still faster than 120,000 lapply calls.
all_reads_dt[, pA_length_stat := ifelse(pA_status == "no_pA", 0, pA_len)]
all_reads_dt[, exonic.length := sapply(strsplit(exon.lengths, ",", fixed =TRUE), function(x) sum(as.numeric(x)))]

dhs.reads.statistics.pA <- all_reads_dt[,
    .(
        ## Calculate summaries for exonic length
        mean.exonic.length = mean(exonic.length),
        median.exonic.length = as.double(median(exonic.length)),
        ## Calculate summaries for exon number
        mean.n.exon = mean(n.exon),
        median.n.exon = as.double(median(n.exon)),
        ## Calculate summary for pA length (using the 0-padded column)
        median.pA.length = as.double(median(pA_length_stat)),
        ## Calculate summaries for genomic length (end - start)
        mean.genomic.length = mean(end - start),
        median.genomic.length = as.double(median(end - start))
    ),
    by = .(geno, pA_status, dhs.type, strand, group_cat) 
]
dhs.reads.statistics.pA <- dhs.reads.statistics.pA[dhs.type %in% dhs.list.filtered, ]

dhs.reads.statistics.merged <- all_reads_dt[,
    .(
        ## Calculate summaries for exonic length
        mean.exonic.length = mean(exonic.length),
        median.exonic.length = as.double(median(exonic.length)),
        ## Calculate summaries for exon number
        mean.n.exon = mean(n.exon),
        median.n.exon = as.double(median(n.exon)),
        ## Calculate summary for pA length (using the 0-padded column)
        median.pA.length = as.double(median(pA_length_stat)),
        ## Calculate summaries for genomic length (end - start)
        mean.genomic.length = mean(end - start),
        median.genomic.length = as.double(median(end - start))
    ),
    by = .(geno,  dhs.type, strand, group_cat) 
]
dhs.reads.statistics.merged <- dhs.reads.statistics.merged[dhs.type %in% dhs.list.filtered, ]

print("Read statistics complete.")

print("Calculating 'nopA' (merged) counts for DESeq2...")

# read counts per replicate, pA+ and pA- together (DESeq2 input)
# It's the same logic as step 3, but without 'pA_status'
dhs.reads.counts.nopA <- all_reads_dt[, .N, by = .(rep, dhs.type, strand, group_cat)]

fwrite(dhs.reads.counts.nopA, file.path(out_dir, "DHS_reads_counts_nopA.tsv"), sep = "\t")


dhs.reads.numbers.nopA <- dhs.reads.counts.nopA
rm(dhs.reads.counts.nopA, dhs.reads.counts.pA, dhs.reads.counts.filtered) # Clean up
print("All processing complete.")

### write important DF
fwrite(dhs.reads.number.summarized.pA, file.path(out_dir, "DHS_reads_number_summarized_pA.tsv"), sep = "\t")
fwrite(dhs.reads.statistics.pA, file.path(out_dir, "DHS_reads_statistics_pA.tsv"), sep = "\t")
fwrite(dhs.reads.statistics.merged, file.path(out_dir, "DHS_reads_statistics_merged.tsv"), sep = "\t")
fwrite(dhs.reads.number.summarized.merged, file.path(out_dir, "DHS_reads_number_summarized_merged.tsv"), sep = "\t")

#### prepare for DESeq2
dhs.prep.FC <- dhs.reads.numbers.nopA %>%
    dplyr::select(dhs.type, group_cat, n = N, rep) %>%
    mutate(dhs.type.group_cat = paste0(dhs.type,".", group_cat)) %>% 
    dplyr::filter(!grepl("ZC3H4", rep)) %>% 
    reshape2::dcast(dhs.type.group_cat~ rep, value.var = "n", fill = 0) %>% 
    column_to_rownames("dhs.type.group_cat") ## %>%
    ## setNames(paste0("scale_",colnames(.)))

df.meta <- colnames(dhs.prep.FC) %>%
    data.frame(sample = .) %>%
    mutate(geno = sample %>% str_replace("mES_","") %>% str_replace("_rep.",""))  %>%
    column_to_rownames("sample") %>%
    mutate(geno = factor(geno, c("WT", "ARS2",  "RBM7",  "ZC3H18", "ZCCHC8")))

deseq.results <- run_deseq(dhs.prep.FC, df.meta)

saveRDS(deseq.results, file.path(out_dir, "deseq_results_DHS.rds"))

## plot log2FC 
deseq.results.df <- deseq.results$df_lfc_shrunk  %>%
    dplyr::select(DHS, log2FoldChange , metric, padj) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    mutate(geno = factor(geno, levels = c("ARS2", "ZC3H18",  "ZCCHC8", "RBM7")))

g <- ggplot(deseq.results.df, aes(x=log2FoldChange, colour = geno )) +
    geom_density() +
    theme_bw() +
    xlim(-1,  2) +
    scale_colour_manual(values = c("#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
    xlab("log2 Fold Change vs WT") +
    ggtitle("DHS log2FC upon exosome related factor KD")
ggsave(g, filename = file.path(out_dir, "DHS_log2FC_density_plot.png"), width = 6, height =4)

### establish list of DHS that are upregulated in each condition
deseq.results.df_up <- deseq.results.df %>%
    drop_na %>% 
    dplyr::filter(log2FoldChange > 1, padj < 0.05 ) %>%
    dplyr::select(DHS, geno) %>%
    mutate(status = "upredulated") %>%
    separate(DHS, into = c("DHS", "type", "group_cat"), sep = "\\.")

## prepare pheatmap

list.dhs <- dhs.reads.number.summarized.merged %>%
    mutate(dhs.type.group_cat = paste(dhs.type, group_cat, sep=".")) %>%
    dplyr::select(dhs.type.group_cat, geno, scaled_total_per_lib, strand, group_cat) %>%
    dplyr::filter(geno != "ZC3H4") %>%
    mutate(DHS = str_replace(dhs.type.group_cat, "\\..*$", "")) %>%
    dplyr::filter(scaled_total_per_lib > 1) %>%
    group_by( DHS) %>%
    add_count() %>%
    ungroup %>%
    dplyr::filter(n > 50) %>%
    .$dhs.type.group_cat%>% unique

deseq.results.df.pheatmap <- deseq.results.df %>%
    dplyr::filter(DHS %in% list.dhs) %>%
    dplyr::select(dhs.type = DHS, log2FoldChange , geno) %>%
    mutate(log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange) %>%
               ifelse(. > 6 , 6, .) %>% ifelse(. < -3, -3, .)) %>%
    separate(dhs.type, into=c("DHS","type","group_cat"), sep="\\.") %>%
    dplyr::filter(type != "gene") %>%
    mutate(geno.type.group_cat= paste(type, geno, group_cat, sep=".")) %>%
    reshape2::dcast(DHS ~ geno.type.group_cat, value.var = "log2FoldChange", fill = 0)  %>%
    column_to_rownames("DHS")

dist.mat <- dist(deseq.results.df.pheatmap)
h <- hclust(dist.mat)

DHS.order <- h$labels[h$order]
col.pal <- c(colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu"))[1:4])(33),
             colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu"))[4:7])(66))

order <- c(2,5,6,7,8,3,4,1) %>%  ##  %>% rep(8) %>% matrix(ncol=8, byrow=TRUE) %>%
    ## t %>%
    ## as.data.frame() %>%
    ## as.list() %>%
    lapply(seq(0, 88, by = 8), function(x, y) { y + x }, .) %>%
    unlist

p <- pheatmap(deseq.results.df.pheatmap[DHS.order, order], 
              cluster_rows = FALSE, 
              cluster_cols = FALSE,
              show_rownames = FALSE,
              breaks = seq(-3,6, length.out = 100),
              color = col.pal,
              width = 6, height = 10)
ggsave(file.path(out_dir,  "sensitivity_heatmap.pdf"), p, width=12, height=8)

###
group_cat.levels <- c("0-100", "100-150", "150-200", "200-250", "250-300", "300-350", "350-400",
                      "400-500", "500-600", "600-700",  "700-800" , "800-1000", "1000-1500",  ">1500")
#### Now we are going to look at genomic length
##############################################

### plot pA length

dhs.pA.length.df <- dhs.reads.statistics.pA %>%
    dplyr::filter(geno != "ZC3H4") %>%
    separate(dhs.type, into = c("DHS", "type"), sep = "\\.") %>%
    dplyr::filter(grepl( "prompt|pmt", type)) %>% 
    mutate(geno = factor(geno, levels = c("WT", "ARS2",  "RBM7",  "ZC3H18", "ZCCHC8"))) %>%
    dplyr::select(DHS, geno, type,  strand, median.pA.length, pA_status, group_cat)  %>%
    dplyr::filter(pA_status != "no_pA") %>%
    mutate(group_cat = factor(group_cat,
                              levels = group_cat.levels)) 

g <- ggplot(dhs.pA.length.df, aes(x = median.pA.length, colour = geno)) +
    geom_density() +
    theme_bw() +
    facet_grid(vars(type), vars(group_cat)) +
    scale_colour_manual(values = c('#828282', "#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
    xlab("Median pA length") +
    ggtitle("DHS read pA length distribution")
ggsave(g, filename = file.path(out_dir, "DHS_read_pA_length_density_plot.png"), width = 10, height =3)

### plot pA fraction
dhs.pA.fraction.df <- dhs.reads.number.summarized.pA %>%
    dplyr::filter(geno != "ZC3H4") %>%
    separate(dhs.type, into = c("DHS", "type"), sep = "\\.") %>%
    dplyr::filter(grepl( "prompt|pmt", type)) %>% 
    mutate(geno = factor(geno, levels = c("WT", "ARS2",  "RBM7",  "ZC3H18", "ZCCHC8"))) %>%
    dplyr::select(DHS, geno, type,   scaled_total_per_lib, pA_status, group_cat)   %>%
    reshape2::dcast(DHS + geno + type + group_cat ~ pA_status, value.var = "scaled_total_per_lib", fill = 0) %>%
    mutate(fraction = pA / (pA + no_pA)) %>%
    mutate(group_cat = factor(group_cat, levels = group_cat.levels)) %>%
    mutate(geno.type  = paste(type, geno, sep = "_"))

for (elt in c("pmt", "prompt")) {
    dhs.pA.fraction.df.elt <- dhs.pA.fraction.df %>%
        dplyr::filter(type == elt)
    g <- ggplot(dhs.pA.fraction.df.elt, aes(x = fraction, colour = geno, fill = geno,
                                        y = after_stat(count))) +
        geom_histogram(binwidth = 0.05) +
        theme_bw() +
        facet_grid(row = vars(geno), cols = vars(group_cat)) +
        scale_colour_manual(values = c('#828282', "#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
        scale_fill_manual(values = c('#828282', "#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
        xlab("pA read fraction") +
        ggtitle(paste0("DHS read pA fraction distribution - ", elt, " reads"))
    ggsave(g, filename = file.path(out_dir, paste0("DHS_read_pA_fraction_histogram_", elt, ".png")), width = 10, height =10)
}

####### Now, we are going to remake some of  these plots but only for the upregulated DHS #######

dhs.up.pA.fraction.df <- dhs.reads.number.summarized.pA %>%
    dplyr::filter(geno != "ZC3H4") %>%
    separate(dhs.type, into = c("DHS", "type"), sep = "\\.") %>%
    dplyr::filter(grepl( "prompt|pmt", type)) %>% 
    mutate(geno = factor(geno, levels = c( "ARS2",  "RBM7",  "ZC3H18", "ZCCHC8"))) %>%
    dplyr::select(DHS, geno, type,   scaled_total_per_lib, pA_status, group_cat)   %>%
    reshape2::dcast(DHS + geno + type + group_cat ~ pA_status, value.var = "scaled_total_per_lib", fill = 0) %>%
    mutate(fraction = pA / (pA + no_pA)) %>%
    mutate(group_cat = factor(group_cat,
                              levels = group_cat.levels)) %>%
    mutate(geno.type  = paste(type, geno, sep = "_")) %>% 
    left_join(deseq.results.df_up, by = c("DHS", "geno", "group_cat", "type")) %>%
    dplyr::filter(!is.na(status)) %>%
    dplyr::filter(status == "upredulated") %>%
    mutate(group_cat =
               factor(group_cat,
                      levels = group_cat.levels)) 

for (elt in c("pmt", "prompt")) {
    dhs.up.pA.fraction.df.elt <- dhs.up.pA.fraction.df %>%
        dplyr::filter(type == elt)
    g <- ggplot(dhs.up.pA.fraction.df.elt, aes(x = fraction, colour = geno, fill = geno,
                                        y = after_stat(count))) +
        geom_histogram(binwidth = 0.05) +
        theme_bw() +
        facet_grid(cols = vars(geno), rows = vars(group_cat)) +
        scale_colour_manual(values = c( "#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
        scale_fill_manual(values = c( "#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
        xlab("pA read fraction") +
        ggtitle(paste0("DHS read pA fraction distribution - ", elt, " reads"))
    ggsave(g, filename = file.path(out_dir, paste0("DHS_up_read_pA_fraction_histogram_", elt, ".png")), width = 6, height =6)
}

#############
#############

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

## deseq.results.df.joined.summarized %>%
##     dplyr::filter(geno =="ZC3H18") %>%
##     dplyr::filter(type == "prompt", group_cat == "800-1000") 

## deseq.results.df.joined %>%
##     dplyr::filter(geno =="ZC3H18") %>%
##     dplyr::filter(type == "prompt", group_cat == "800-1000", pA_status == "pA")  %>%
##     head

deseq.results.df.joined.summarized <- deseq.results.df.joined %>%
    ## group_by(geno, group_cat, pA_status) %>%
    ## add_count() %>% 
    group_by(geno, group_cat, status, pA_status, type) %>%
    dplyr::summarize(  n = n() ) %>%
    group_by(geno, group_cat, pA_status, type) %>%
    mutate( n.total = sum(n), fraction = n/n.total) %>%
    mutate(group_cat =
               factor(group_cat,
                      levels = group_cat.levels)) %>%
    mutate(geno = factor(geno, levels = c("ARS2", "RBM7", "ZC3H18",  "ZCCHC8")))  

for (elt in c("pmt", "prompt", "gene", "exon1")) {
    deseq.results.df.joined.summarized.elt <- deseq.results.df.joined.summarized %>%
        dplyr::filter(type == elt)
    g <- ggplot(deseq.results.df.joined.summarized.elt, aes(x = geno, y = n, fill = geno, alpha = status)) +
        geom_bar(stat ="identity" , color = "black", size = 0.2) + #, position = "dodge") +
        ## scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x  = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 6),
              strip.placement = "outside") +
        facet_grid(vars( pA_status), vars(group_cat), switch = "x") +
        scale_fill_manual(values = c("#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
        scale_alpha_manual(values = c("upregulated" = 1, "unchanged" = 0.6)) +
        ## scale_color_manual(values = c("upregulated" = "white", "unchanged" = "grey")) +
        ## scale_fill_manual(values = c("upregulated" = "#E07A5F", "total" = "grey")) +
        xlab("DHS length category") +
        ylab("Count of DHS") +
        ggtitle(paste0("DHS upregulated count upon exosome related factor KD - ", elt, " reads"))
    ggsave(g, filename = file.path(out_dir, paste0("DHS_up_count_barplot_length_separated_pA_nopA_", elt, ".png")), width = 10, height =4)
    g <- ggplot(deseq.results.df.joined.summarized.elt %>%
              dplyr::filter(status =="upregulated"), aes(x = geno, y = fraction, fill = geno)) +
        geom_bar(stat ="identity" , color = "black", size = 0.2) + #, position = "dodge") +
        ## scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x  = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 6),
              strip.placement = "outside") +
        facet_grid(vars( pA_status), vars(group_cat), switch = "x") +
        scale_fill_manual(values = c("#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
        ## scale_color_manual(values = c("upregulated" = "white", "unchanged" = "grey")) +
        ## scale_fill_manual(values = c("upregulated" = "#E07A5F", "total" = "grey")) +
        xlab("DHS length category") +
        ylab("fraction of upregulated DHS") +
        ggtitle(paste0("DHS upregulated count upon exosome related factor KD - ", elt, " reads"))
    ggsave(g, filename = file.path(out_dir, paste0("DHS_up_fraction_barplot_length_separated_pA_nopA_", elt, ".png")), width = 10, height =4)
}

## plot logFC in violin plots
deseq.results.df.joined.plot.violin <-  deseq.results.df.joined %>%
    dplyr::select(-status, -tpm.pA, -tpm.merged, -pA_status) %>%
    distinct() %>%
    mutate(group_cat =
               factor(group_cat,
                      levels = group_cat.levels)) %>%
    mutate(geno = factor(geno, levels = c("ARS2", "RBM7", "ZC3H18",  "ZCCHC8")))

for (elt in c("pmt", "prompt", "gene", "exon1")) {
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
        scale_fill_manual(values = c("#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
        scale_alpha_manual(values = c("upregulated" = 1, "unchanged" = 0.6)) +
        xlab("DHS length category") +
        ylab("log2 Fold Change vs WT") +
        ggtitle(paste0("DHS log2FC upon exosome related factor KD - ", elt, " reads"))
    ggsave(g, filename = file.path(out_dir, paste0("DHS_log2FC_violin_plot_length_separated_", elt, ".png")), width = 10, height =4)
}

g <- ggplot(deseq.results.df.joined.plot.violin, aes( x = group_cat, y= log2FoldChange, fill = geno)) +
    geom_violin( color = "black", size = 0.2) +
    theme_bw() +
    theme(# axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x  = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside") +
    stat_summary(fun =mean, geom="point", shape=1, size=2, color="black") +
    stat_summary(fun = median, geom="point", shape=2, size=2, color="black") +
    facet_grid(row = vars(geno),  switch = "x") +
    scale_fill_manual(values = c("#E07A5F", "#3D405B", "#81B29A", "#F2CC8F")) +
    scale_alpha_manual(values = c("upregulated" = 1, "unchanged" = 0.6)) +
    xlab("DHS length category") +
    ylab("log2 Fold Change vs WT") +
    ggtitle("DHS log2FC upon exosome related factor KD")
ggsave(g, filename = file.path(out_dir, "DHS_log2FC_violin_plot_length_separated.png"), width = 10, height =4)


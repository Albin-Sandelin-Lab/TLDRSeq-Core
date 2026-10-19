## pA+ and pA- PROMPTs in WT, RBM7 and ZCCHC8: 3' end peaks, sequence logos, motifs
## (XSTREME), motif positions, 3' end shift between pA+ and pA- reads.
##
## Usage (from 05-PROMPT/): Rscript 05_prompt_pA_motifs.r
##
## PROMPT reads here: reads of the DHS on the strand opposite to the gene that do not
## overlap the gene body (GENCODE vM35 gene, any strand).
## A DHS has a pA+ (pA-) PROMPT in a genotype if it has > 5 pA+ (pA-) PROMPT reads.
## Clusters: prompts.WT (pA+ PROMPT in WT), prompts.NoWT.ZCCHC8 / prompts.NoWT.RBM7
## (pA+ PROMPT in the KD but not in WT).
##
## Output: $PROMPT_DIR/prompt_pA_motifs/
##   pA vs total read count plots, Venn / UpSet of PROMPT sets
##   seq/, seq_for_motifs/   sequences around PROMPT 3' end peaks (FASTA)
##   logo2_all_pA{+,-}.pdf, nt_content.pdf
##   xstreme2/<cluster>/     XSTREME motifs and PFMs
##   Freq_position_PWM_per_cluster.pdf  position of the best motif hits
##   median_diff_*.pdf       3' end shift between pA+ and pA- reads of the same PROMPT
##
## PROMPT reads overlapping the gene body are excluded here (unlike 01), on purpose.
##
## Requires MEME suite (xstreme) and lib/meme2pfm2.sh.
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(parallel)
library(pheatmap)
library(data.table)
library(plyranges)
library(org.Mm.eg.db)
library(tidyverse)
library(ggVennDiagram)

source("../config/paths.R")
source("../lib/tldr_reads.R")
source("../lib/peak_calling.r")
source("../lib/pwm_functions.r")

genome <- Biostrings::readDNAStringSet(GENOME_FA, format = "fasta", use.names = TRUE)

gtf <- import(GENCODE_GTF, format = "GTF")
gtf.reduced <-  gtf %>%
    plyranges::select(type, gene_id, gene_type)
    

########################
##### new file to compare
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

df_DHS_prompt_map <- df_DHS_prompt_map %>%
    dplyr::select(DHS) %>%
    left_join(dhs.gene.strand) %>%
    dplyr::select(DHS, gene_id =gencode.intersect, gene_strand) %>%
    mutate(direction_prompt = ifelse(gene_strand == "+" , "-", "+"), direction_DHS = gene_strand)

#################### Code
out_dir <- file.path(PROMPT_DIR, "prompt_pA_motifs")
dir.create(out_dir, showWarnings=FALSE, recursive = TRUE)

genos <- c("WT", "ZCCHC8", "RBM7")

## Reads of the PROMPT of one DHS in one genotype.
## extraction.mode = TRUE : returns c(n pA+ reads, n pA- reads, strand, DHS, genotype)
## extraction.mode = FALSE: returns list(pA+ / pA- read counts per 3' end position, strand)
get.signal.as.df2<- function(dhs, geno, extraction.mode = FALSE) {
    ## print(dhs)
    dhs_info <- df_DHS_prompt_map %>% 
        filter(DHS == dhs)
    dir_prompt <- dhs_info$direction_prompt %>% as.character()
    ##### this has been replaced by the line below
    ## gene_id.entrez <- dhs_info$gene_id %>% as.character() ## 
    ## gene_id.ensembl <- ensembl_ids[names(ensembl_ids) == gene_id.entrez]  %>% str_replace("\\..*","")
    #######
    gene_id.ensembl <- dhs_info$gene_id 
    if(is.na(gene_id.ensembl)) {
        if(extraction.mode) {
            extraction.c <- c(0, 0, dir_prompt, dhs, geno)
            return(extraction.c)
        }
        return(NULL)
    }
    ## get reads
    tmp2 <- tmp_replicates[grepl(geno, tmp_replicates)] 
    DHS_reads = read_TLDR_reads_by_list_of_DHS(query_list_DHS = dhs, query_list_replicates = tmp2)
    if(nrow(DHS_reads) == 0) {
        if(extraction.mode) {
            extraction.c <- c(0, 0, dir_prompt, dhs, geno)
            return(extraction.c)
        }
        return(NULL)
    }
    DHS_reads <- DHS_reads %>%
        mutate(genotype = str_replace(V14, ".*_(.*)_.*","\\1")) %>%
        dplyr::filter(geno == genotype) %>%    
        dplyr::select(seqnames = "V1", start = "V2", end = "V3", pA_positive = "V5", strand = "V6")  %>%
        ## dplyr::filter(strand == dir_prompt)
        as_granges %>%
        .[strand(.) == dir_prompt]
    gene.location <- gtf.reduced[gtf.reduced$type =="gene"] %>%
        .[grepl(gene_id.ensembl, .$gene_id)] 
    DHS_reads.out <- subsetByOverlaps(DHS_reads, gene.location, invert = TRUE, ignore.strand = TRUE)
    ## n.overlap <- subsetByOverlaps(DHS_reads, gene.location, ignore.strand = TRUE)
    ## DHS_reads.out <- DHS_reads ## debug
    if (extraction.mode) {
        extraction.c <- c(length(DHS_reads.out[DHS_reads.out$pA_positive > 15]),
                          length(DHS_reads.out[DHS_reads.out$pA_positive <= 15]),
                          dir_prompt,
                          dhs,
                          geno)
        return(extraction.c)
    }##
    DHS.3prime <- DHS_reads.out  %>%
        anchor_3p() %>%
        mutate(width = 1)  %>%
        as.data.frame %>%
        dplyr::select(TTS = start, pA_positive) %>% 
        mutate(pA_positive = factor(ifelse(pA_positive <= 15 , "nopA", "pA"), levels = c("pA", "nopA"))) %>%
        group_by(TTS, pA_positive) %>%
        dplyr::summarize(n = n(), .groups = "drop" ) %>%
        reshape2::dcast(TTS ~ pA_positive, value.var = "n", fill = 0, drop = FALSE)
    return(list(DHS.3prime, dir_prompt))
}

## 3' end peaks of the pA+ and pA- PROMPT reads of one DHS (find_peaks, lib/peak_calling.r)
get.peaks<- function(dhs, geno) {
    ## print(dhs)
    dhs_info <- df_DHS_prompt_map %>% 
        filter(DHS == dhs)
    dir_prompt <- dhs_info$direction_prompt %>% as.character()
    ## this has been replaced by the line below
    ## gene_id.entrez <- dhs_info$gene_id %>% as.character()
    ## gene_id.ensembl <- ensembl_ids[names(ensembl_ids) == gene_id.entrez]  %>% str_replace("\\..*","")
    #######
    gene_id.ensembl <- dhs_info$gene_id 
    if (is.na(gene_id.ensembl)) {
        return(NULL)
    }
    ## get reads
    tmp2 <- tmp_replicates[grepl(geno, tmp_replicates)] 
    DHS_reads = read_TLDR_reads_by_list_of_DHS(query_list_DHS = dhs, query_list_replicates = tmp2)
    if(nrow(DHS_reads) == 0) {
        return(NULL)
    }
    DHS_reads <- DHS_reads %>%
        dplyr::filter(V6 == dir_prompt) %>%
        dplyr::filter(grepl(geno,  V14)) %>%    
        dplyr::filter(V6 == dir_prompt)
    gene.location <- gtf.reduced[gtf.reduced$type =="gene"] %>%
        .[grepl(gene_id.ensembl, .$gene_id)]
    DHS_reads.gr <- GRanges(DHS_reads$V1,  IRanges(DHS_reads$V2, DHS_reads$V3),    strand = DHS_reads$V6)
    DHS_reads.ov <- findOverlaps(DHS_reads.gr, gene.location, ignore.strand = TRUE)
    DHS_reads.ov.ind <- queryHits(DHS_reads.ov) %>% unique
    DHS_reads.out.ind <- setdiff(1:nrow(DHS_reads), DHS_reads.ov.ind)
    DHS_reads.out <- DHS_reads[DHS_reads.out.ind,]
    if (nrow(DHS_reads.out) == 0) {
        return(NULL)
    }
    if(sum(DHS_reads.out$V5 > 15) != 0) {
        peaks.pA <- find_peaks(DHS_reads.out, strand = dir_prompt,window_count_filter = -1, kde_bw = 5,
                               kde_min_gradient_length=2,
                               pA="+")
    } else {
        peaks.pA <-  data.frame()
    }
    if(sum(DHS_reads.out$V5 <= 15) != 0) {
        peaks.nopA <- find_peaks(DHS_reads.out, strand = dir_prompt,window_count_filter = -1, kde_bw = 5,
                               kde_min_gradient_length=2,
                               pA="-")
    } else {
        peaks.nopA <- data.frame()
    }
    peaks <- rbind.data.frame(peaks.pA, peaks.nopA) %>%
        mutate(dhs = dhs)
    if (nrow(peaks) == 0) {
        return(NULL)
    }
    return(peaks)
}


extraction.list <- list()
for (elt in genos)
{
    extraction.list[[elt]] <- mclapply(X = df_DHS_prompt_map$DHS,
                                       FUN = function(x, y) {
                                           get.signal.as.df2(dhs = x, geno = y, extraction.mode = TRUE)},
                                       elt,
                                       mc.cores = 70)
}

extraction.df <- list()
for (elt in genos)
{
    extraction.df[[elt]] <- unlist(extraction.list[[elt]]) %>%
        {data.frame(pA = .[rep(c(TRUE, FALSE, FALSE, FALSE, FALSE), length(.)/5)],
                    nopA = .[rep(c(FALSE, TRUE, FALSE, FALSE,  FALSE), length(.)/5)],
                    strand = .[rep(c(FALSE, FALSE, TRUE,FALSE, FALSE), length(.)/5)],
                    name = .[rep(c(FALSE, FALSE, FALSE, TRUE, FALSE), length(.)/5)],
                    geno = elt)}
}

    

extraction.df <- do.call(rbind.data.frame, extraction.df)  %>%
    mutate(pA = as.numeric(pA),
           nopA = as.numeric(nopA))

df.dhs <- extraction.df %>%
    split(f=.$geno) %>%
    lapply(remove_rownames) %>%
    lapply(column_to_rownames, "name") %>%
    lapply(function(x) setNames(x, paste0(colnames(x),".", x$geno[1]))) %>%
    lapply(rownames_to_column, "name") %>% 
    {left_join(left_join(.[[1]], .[[2]], by = "name"), .[[3]], by = "name")} %>%
    dplyr::select(name, strand = strand.ZCCHC8, pA.RBM7, nopA.RBM7, pA.WT, nopA.WT, pA.ZCCHC8, nopA.ZCCHC8) 

#### make x y plots
g <- ggplot(df.dhs %>% dplyr::filter(pA.WT < 5 ), aes(y = pA.RBM7, x = pA.RBM7 + nopA.RBM7)) +
    geom_hex() +
    scale_fill_viridis_c(trans = "log10")+
    ggpubr::stat_cor(method = "pearson") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "pA+ Reads Counts RBM7 vs all counts RBM7 (no WT pA+ prompts in these)",
         x = "Total read  Counts RBM7",
         y = "pA+ Reads Counts RBM7") +
    theme_bw()
ggsave(file.path(out_dir, "pA_RBM7_vs_total_RBM7_WT_filt.pdf"), g, width = 8, height = 6)

g <- ggplot(df.dhs %>% dplyr::filter(pA.WT < 5 ), aes(y = pA.ZCCHC8, x = pA.ZCCHC8 + nopA.ZCCHC8)) +
    geom_hex() +
    scale_fill_viridis_c(trans = "log10")+
    ggpubr::stat_cor(method = "pearson") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "pA+ Reads Counts ZCCHC8 vs all counts ZCCHC8 (no WT pA+ prompts in these)",
         x = "Total read  Counts ZCCHC8",
         y = "pA+ Reads Counts ZCCHC8") +
    theme_bw()
ggsave(file.path(out_dir, "pA_ZCCHC8_vs_total_ZCCHC8_WT_filt.pdf"), g, width = 8, height = 6)

custom_blues <- RColorBrewer::brewer.pal(9, "Blues")[3:9]
g <- ggplot(df.dhs, aes(y = pA.WT, x = pA.WT + nopA.WT)) +
    geom_hex() +
    ## scale_fill_viridis_c(trans = "log10")+
    scale_fill_gradientn(colors = custom_blues, trans = "log10") +
    ggpubr::stat_cor(method = "pearson") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "pA+ Reads Counts WT vs WT",
         x = "Total read  Counts WT",
         y = "pA+ Reads Counts WT") +
    theme_bw()
ggsave(file.path(out_dir, "pA_WT_vs_total_WT.pdf"), g, width = 8, height = 6)

##### We will remake the same plots, but we are going to change the color scales. (V2, using RdYlBu)

custom_blues <- RColorBrewer::brewer.pal(9, "Blues")[3:9]
g <- ggplot(df.dhs %>% dplyr::filter(pA.WT < 5 ), aes(y = pA.RBM7, x = pA.RBM7 + nopA.RBM7)) +
    geom_hex() +
    ## Feed the truncated palette into scale_fill_gradientn
    scale_fill_gradientn(colors = custom_blues, trans = "log10") +
    ggpubr::stat_cor(method = "pearson") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "pA+ Reads Counts RBM7 vs all counts RBM7 (no WT pA+ prompts in these)",
         x = "Total read Counts RBM7",
         y = "pA+ Reads Counts RBM7") +
    theme_bw()
## Updated filename
ggsave(file.path(out_dir, "pA_RBM7_vs_total_RBM7_WT_filt_Blues.pdf"), g, width = 8, height = 6)

g <- ggplot(df.dhs %>% dplyr::filter(pA.WT < 5 ), aes(y = pA.ZCCHC8, x = pA.ZCCHC8 + nopA.ZCCHC8)) +
    geom_hex() +
    ## Feed the truncated palette into scale_fill_gradientn
    scale_fill_gradientn(colors = custom_blues, trans = "log10") +
    ggpubr::stat_cor(method = "pearson") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "pA+ Reads Counts ZCCHC8 vs all counts ZCCHC8 (no WT pA+ prompts in these)",
         x = "Total read Counts ZCCHC8",
         y = "pA+ Reads Counts ZCCHC8") +
    theme_bw()
## Updated filename
ggsave(file.path(out_dir, "pA_ZCCHC8_vs_total_ZCCHC8_WT_filt_Blues.pdf"), g, width = 8, height = 6)

### for Venn diagram
prompts.ZCCHC8  <- df.dhs   %>% dplyr::filter(pA.ZCCHC8 > 5) %>% .$name
prompts.WT  <- df.dhs   %>% dplyr::filter(pA.WT > 5) %>% .$name
prompts.RBM7 <- df.dhs   %>% dplyr::filter(pA.RBM7 > 5) %>% .$name
prompts.nopA.ZCCHC8  <- df.dhs   %>% dplyr::filter(nopA.ZCCHC8 > 5) %>% .$name
prompts.nopA.WT  <- df.dhs   %>% dplyr::filter(nopA.WT > 5) %>% .$name
prompts.nopA.RBM7 <- df.dhs   %>% dplyr::filter(nopA.RBM7 > 5) %>% .$name

x <- list("ZCCHC8 pA+" = prompts.ZCCHC8,
          "WT pA+" = prompts.WT,
          "ZCCHC8 pA-" = prompts.nopA.ZCCHC8,
          "RBM7 pA+" = prompts.RBM7,
          "RBM7 pA-" = prompts.nopA.RBM7)

g <- ggVennDiagram::ggVennDiagram(x, 
                       label_alpha = 0.5,
                       label = "count",
                       set_size = 3) +
    ggplot2::ggtitle("Prompts") +
    scale_fill_gradient(low="grey90",high = "red") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    scale_x_continuous(expand = expansion(mult = .05))
ggsave(file.path(out_dir, "venn_diagram_ggplot.pdf"), g, width = 7, height = 7) 

y <- x
y[[6]] <- prompts.nopA.WT
names(y) <- c(names(x), "WT pA-") 
g <- ggVennDiagram::ggVennDiagram(y, force_upset = TRUE, specific =TRUE)
ggsave(file.path(out_dir, "upset_ggplot.pdf"), g, width = 10, height = 4) 

g <- ggVennDiagram::ggVennDiagram(y , 
                                  label_alpha = 0.5,
                                  label = "count",
                                  set_size = 3) +
    ggplot2::ggtitle("Prompts") +
    scale_fill_gradient(low="grey90",high = "red") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    scale_x_continuous(expand = expansion(mult = .05))
ggsave(file.path(out_dir, "venn_diagram_ggplot_all.pdf"), g, width = 7, height = 7) 

#######

#######

dhs.ZCCHC8.only <-  prompts.ZCCHC8[!prompts.ZCCHC8 %in% prompts.WT]
dhs.RBM7.only <-  prompts.RBM7[!prompts.RBM7 %in% prompts.WT]
dhs.WT <- prompts.WT

peaks.ZCCHC8 <-  mclapply(dhs.ZCCHC8.only, get.peaks, "ZCCHC8" , mc.cores = 50)
peaks.ZCCHC8.df <-   do.call(rbind.data.frame, peaks.ZCCHC8)
peaks.WT <-  mclapply(dhs.WT, get.peaks, "WT" , mc.cores = 50)
peaks.WT.df <-   do.call(rbind.data.frame, peaks.WT)
peaks.RBM7 <-  mclapply(dhs.RBM7.only, get.peaks, "RBM7" , mc.cores = 50)
peaks.RBM7.df <-   do.call(rbind.data.frame, peaks.RBM7)

list.peaks.df <- list(prompts.NoWT.ZCCHC8 = peaks.ZCCHC8.df, 
                      prompts.WT = peaks.WT.df,
                      prompts.NoWT.RBM7 = peaks.RBM7.df) %>%
    mapply(mutate, ., cluster = names(.), SIMPLIFY = FALSE)  %>%
    do.call(rbind.data.frame,.)

format.peaks <- function(dhs.peaks, th = 10, max.peak = FALSE)
{
    maximas.values <- str_split(dhs.peaks$local_maxima_window_counts, pattern = ",") %>%
        lapply(as.numeric)
    global.maxima.index <- lapply(maximas.values, which.max)
    peaks <- str_split(dhs.peaks$local_maxima, pattern = ",") %>%
        lapply(as.numeric)
    peaks.adj <- str_split(dhs.peaks$local_maxima.adj, pattern = ",") %>%
        lapply(as.numeric)
    strand <- dhs.peaks$strand
    pA_status <- dhs.peaks$pA
    cluster <- dhs.peaks$cluster
    chr <- str_replace(dhs.peaks$dhs, "_.*", "")
    if(max.peak)
    {
    filtered_peaks <- mapply(function(x, y, z) x[z][y[z] >= th], peaks.adj, maximas.values, global.maxima.index)  %>% 
        lapply(function(x) if(length(x) == 0) NA else x)
    }else{
    filtered_peaks <- mapply(function(x, y) x[y >= th], peaks.adj, maximas.values) %>% ## changed to max.adj
        lapply(function(x) if(length(x) == 0) NA else x)
    }
    peaks.bed <- mapply(data.frame, 
                        seqnames = chr, 
                        start = filtered_peaks , 
                        end = filtered_peaks, 
                        strand = strand,
                        pA_status = pA_status,
                        cluster = cluster,
                        SIMPLIFY=FALSE)
    return(peaks.bed)
}

## run for all the peaks at once

#### new addition ######################################
#### we don't want pA+ peaks that overlap with pA- peaks


peaks.all.curated <- format.peaks(list.peaks.df, 5) %>%
    do.call(rbind.data.frame, .) %>%
    drop_na %>%
    plyranges::as_granges() %>%
    plyranges::anchor_3p()  %>%
    mutate(width =  width + 50) %>%
    plyranges::anchor_5p() %>%
    mutate(width = width + 50)  %>% 
    mutate(name = paste0(.$cluster,".", .$pA_status)) %>% 
    split(f= .$name )

peaks.pAplus.curated <- peaks.all.curated[grepl("\\+$", names(peaks.all.curated))]
peaks.pAminus.curated <- peaks.all.curated[grepl("-$", names(peaks.all.curated))]

peaks.pAplus.curated.subset <- mapply(function(x, y) subsetByOverlaps(x, y, invert=TRUE, ignore.strand=TRUE),
                                     peaks.pAplus.curated,
                                     peaks.pAminus.curated,
                                     SIMPLIFY=FALSE) %>%
    as("GRangesList") 

peaks.pAminus.curated.subset <- mapply(function(x, y) subsetByOverlaps(x, y, invert=TRUE, ignore.strand=TRUE),
                                       peaks.pAminus.curated,
                                       peaks.pAplus.curated,                                       
                                       SIMPLIFY=FALSE) %>%
    as("GRangesList") 

peaks.all.logo <- c(peaks.pAplus.curated.subset, peaks.pAminus.curated.subset) %>%
    unlist() %>%
    anchor_3p() %>%
    mutate(width = width + 100) %>%
    subsetByOverlaps(gtf.reduced, invert = TRUE, ignore.strand = TRUE) %>%
    split(f=.$name)

#######
sequences.list <- lapply( peaks.all.logo, function(x, y) BSgenome::getSeq(y, x), genome  ) %>%
    lapply(function(x) x[!grepl("[^ACTG]", as.character(x))]) 

seq_dir <- file.path(out_dir, "seq")
dir.create(seq_dir, showWarnings= FALSE, recursive = TRUE)
a <- sequences.list %>%
    mapply(function(x,y) setNames(x, paste0(y, "_", c(1: length(x)))),. ,  names(.), SIMPLIFY=FALSE) %>% 
    mapply(Biostrings::writeXStringSet, .,  file.path(seq_dir, paste0(names(.), ".fasta")), SIMPLIFY=FALSE)

################ make motifs for each cluster
############################################

###### Now, pot logos of 5' of each cluster
nt.freq <- function(seq.logo, reduce.left = 0, reduce.right = 0)
{
    freq <- seq.logo %>%
        as.character()  %>%
        str_split(., "") %>%
        do.call(rbind, .) %>%
        as.data.frame %>%
        as.list()  %>%
        lapply(paste, collapse = "") %>%
        lapply(function(x) sapply(list(A="A",C="C",G="G",T="T"), function(y) str_count( x ,y ))) %>%
        do.call(cbind.data.frame, .) %>%
        t %>%
        as.data.frame %>%
        remove_rownames %>%
        mutate(position = 1:nrow(.)) %>%
        dplyr::filter(position %in% seq(reduce.left, max(position ) - reduce.right)) %>%
        mutate(position = position - reduce.left) %>% 
        reshape2::melt(id.vars=c("position"))
    return(freq)
}

freq.logo <- mclapply(sequences.list, nt.freq, 100, 0,  mc.cores=4) ## %>%
    ## .[c("prompts.WT.+","prompts.NoWT.ZCCHC8.+", "prompts.NoWT.RBM7.+")] 

freq.nt <- mclapply(sequences.list, nt.freq, mc.cores=4)

nt.content <- sapply(lapply(sequences.list, Biostrings::letterFrequency, "GC"), mean) /200

nt.content.list <- freq.nt %>%
                                        #[c("prompts.WT.+","prompts.NoWT.ZCCHC8.+", "prompts.NoWT.RBM7.+")] %>%
    lapply(function(x) group_by(x, position) %>% mutate(freq = value / sum(value)))  %>%
    mapply(mutate, ., seq_group = names(.), SIMPLIFY = FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames %>%
    as.data.frame %>% 
    mutate(position = position - 150)

g <- ggplot(nt.content.list, aes(x = position, y = freq, color = variable)) +
    geom_line() +
    facet_wrap(~seq_group, nrow = 3) +
    theme_bw()
ggsave(file.path(out_dir, "nt_content.pdf"), g, width = 10, height = 5)
    

all.freq <- mapply(mutate, freq.logo,
                   seq_group = names(freq.logo), SIMPLIFY = FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    dplyr::rename(nt=variable) %>%
    group_by(seq_group, position) %>%
    mutate(n.seq = sum(value)) %>% 
    group_by(nt, seq_group) %>%
    mutate(nt.abundance = sum(value)) %>%
    ungroup %>%
    group_by(seq_group) %>%
    mutate(nt.tot=sum(value)) %>%
    ungroup %>% 
    mutate(freq=value/n.seq, freq.norm = value/(n.seq * (nt.abundance/nt.tot))) 

all.logo2.towrite <- all.freq %>% 
    mutate(freq = freq)  %>%
    dplyr::select(seq_group, position, nt, freq) %>%
    group_by(position, seq_group) %>%
    mutate(freq = freq + 0.01) %>%
    mutate(tot = sum(freq), freq = (freq/tot) ) %>%
    ungroup %>% 
    reshape2::dcast(seq_group  + position ~ nt, value.var = "freq") %>% 
    remove_rownames %>% 
    split(f=.$seq_group) %>%
    lapply(dplyr::select, A,C,G,T)  %>%
    lapply(remove_rownames)

p_list2 <- all.logo2.towrite %>%
    lapply(function(x) t(as.matrix(x) *  matrix(rep(2 + rowSums(( as.matrix(x)) *log (as.matrix(x), 2)),4),ncol=4)))  

g <- ggseqlogo::ggseqlogo(p_list2[grepl("\\+", names(p_list2))] %>% setNames( names(.)),
                          method="custom", seq_type="dna") +
    facet_wrap(~seq_group, ncol=1, scales='free', dir ="v") +
    scale_x_continuous(breaks=seq(1,101, 8), labels=seq(-50, 50 , 8), guide = guide_axis(n.dodge = 2))  
ggpubr::ggexport(g, filename = file.path(out_dir, "logo2_all_pA+.pdf"), width = 10, height = 5 )

g <- ggseqlogo::ggseqlogo(p_list2[grepl("-", names(p_list2))] %>% setNames( names(.)),
                          method="custom", seq_type="dna") +
    facet_wrap(~seq_group, ncol=1, scales='free', dir ="v") +
    ylim(c(0, 0.6)) +
    scale_x_continuous(breaks=seq(1,101, 8), labels=seq(-50, 50 , 8), guide = guide_axis(n.dodge = 2))  
ggpubr::ggexport(g, filename = file.path(out_dir, "logo2_all_pA-.pdf"), width = 10, height = 5 )

###### get prompt contexts and find some motifs.
### get best peak per sequence

peaks.all.max <- format.peaks(list.peaks.df, 5, TRUE) %>%
    do.call(rbind.data.frame, .) %>%
    drop_na %>%
    plyranges::as_granges() %>%
    plyranges::anchor_3p()  %>%
    mutate(width = width + 150) %>%
    plyranges::anchor_5p() %>%
    mutate(width = width  + 50)  %>% 
    subsetByOverlaps(gtf.reduced, invert = TRUE, ignore.strand = TRUE) %>% 
    split(f=paste0(.$cluster,".", .$pA_status))   

### we are going to overlap with the previous filtered peaks to avoid pA+ peaks that overlap with pA- peaks
peaks.pAplus.max <- peaks.all.max[grepl("\\+$", names(peaks.all.max))] 

peaks.pAplus.max.curated <- mapply(function(x, y) subsetByOverlaps(x, y, ignore.strand=FALSE, minoverlap=200),
                                     peaks.pAplus.max,
                                     peaks.pAplus.curated %>% resize(width = 201, fix = "end"),
                                     SIMPLIFY=FALSE) %>%
    as("GRangesList")

#######
sequences.list.3 <- lapply( peaks.pAplus.max.curated %>% setNames( names(.) %>% str_replace("..$", "")),
                           function(x, y) BSgenome::getSeq(y, x), genome  ) %>%
    lapply(function(x) x[!grepl("[^ACTG]", as.character(x))]) 
names.for.latter <- names(sequences.list.3)

seq_dir3 <- file.path(out_dir, "seq_for_motifs")
dir.create(seq_dir3, showWarnings= FALSE, recursive = TRUE)
a <- sequences.list.3 %>%
    mapply(function(x,y) setNames(x, paste0(y, "_", c(1: length(x)))),. ,  names(.), SIMPLIFY=FALSE) %>% 
    mapply(Biostrings::writeXStringSet, .,  file.path(seq_dir3, paste0(names(.), ".fasta")), SIMPLIFY=FALSE)

#### re-read the sequences so that they ghave a name ....
sequences.list.3 <- file.path(seq_dir3, paste0( names(sequences.list.3), ".fasta")) %>%
    lapply(Biostrings::readDNAStringSet, format = "fasta", use.names = TRUE)
names(sequences.list.3) <- names.for.latter

#### all of the comming block of code is there to generate logo
##################################################
bg.sequences.list.3 <- lapply(names(sequences.list.3), function(x) {
  sequences.list.3[names(sequences.list.3) != x] %>% setNames(NULL) %>% do.call(c, .)
}) %>% setNames(names(sequences.list.3))

names.for.latter <- names(sequences.list.3)

a <- bg.sequences.list.3 %>%
    mapply(function(x,y) setNames(x, paste0(y, "_", c(1: length(x)))),. ,  names(.), SIMPLIFY=FALSE) %>% 
    mapply(Biostrings::writeXStringSet, .,  file.path(seq_dir3, paste0("bg_", names(.), ".fasta")), SIMPLIFY=FALSE)

bg.sequences.list <- file.path(seq_dir3, paste0( "bg_", names(sequences.list.3), ".fasta")) %>%
    lapply(Biostrings::readDNAStringSet, format = "fasta", use.names = TRUE)

## de novo motifs for each cluster (pA+ PROMPT 3' ends)
command <- paste("module load meme; xstreme --p",
                 file.path(seq_dir3, paste0( names(sequences.list.3), ".fasta")),
                 "--oc", file.path(out_dir, "xstreme2", paste0(names(sequences.list.3))),
                 "--meme-nmotifs 3 --streme-nmotifs 5 --rna" )
commands <- paste(paste(command, collapse = " & "), "& wait")
system(commands)

## Convert the XSTREME motifs to PFM files (one <motif>.pfm per motif, in the same folder)
meme2pfm <- normalizePath("../lib/meme2pfm2.sh")
list.commands <- c(file.path(out_dir, "xstreme2",  names(sequences.list.3))) %>%
    lapply(list.files, pattern = ".*xstreme.txt$", full.names=TRUE) %>%
    .[sapply(., function(x) length(x) != 0)] %>% 
    lapply(function(x, y) paste(y, x), meme2pfm)

commands <- paste(list.commands, collapse = "; ")
system(commands)

#################
##########

#### get pwm from RBM7 dir
list.pfm.name <- file.path(out_dir, "xstreme2",  "prompts.NoWT.RBM7") %>% 
    lapply(list.files, pattern = ".*pfm$", full.names=TRUE) 

pfm.names <-    file.path(out_dir, "xstreme2",  "prompts.NoWT.RBM7") %>% 
    lapply(list.files, pattern = ".*pfm$", full.names=TRUE)  %>%
    lapply(basename) %>%
    lapply(str_replace, pattern = "\\.pfm$", replacement = "")


list.pfm <- lapply(list.pfm.name, function(x) lapply( x, read.table, header=TRUE, sep="\t", skip=1)) 

list.pwm <- lapply(list.pfm, function(x) lapply(x, make.PWM))
list.pwm <- mapply(function(x, y, zz) {
    z <- setNames(x, paste0(y))
    return(z)},
    list.pwm, pfm.names, SIMPLIFY=FALSE)

list.pwm.unlisted <- list(unlist(list.pwm, recursive = FALSE))

#### run PWM on sequences
#########################
sequences.list.char <- sequences.list.3 %>%
    mclapply(as.character, mc.cores=12)

##### now, locate the best score of each sequence

sequence.best_pos.scores <- mcmapply(run_all_pwm_on_seq, rep(list.pwm.unlisted, length(sequences.list.char)),
                                     sequences.list.char,
                                     TRUE,
                                     0.2,
                                     SIMPLIFY=FALSE, mc.cores=2)
names(sequence.best_pos.scores) <- names(sequences.list.char)

position <- transpose(sequence.best_pos.scores)

position.df <- lapply(position, function(x) {
    mapply(function(y, z, zz) as.data.frame(table( y = c(y, 1:zz)))  %>% mutate(Freq = (Freq -1)/z),
           x, lapply(sequences.list.char, length), MoreArgs = list(zz = nchar(sequences.list.char[[1]][[1]])), 
           SIMPLIFY=FALSE)}) %>%
    lapply(., function(x) mapply( function(y, name) mutate(y, cluster = name), x, names(x) ,  SIMPLIFY=FALSE)) %>%
    lapply(function(x) do.call( rbind.data.frame, x)) %>%
    mapply(mutate, ., PWM = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    mutate(y = y %>% as.character %>% as.numeric)
    

position.matrix.list <- split(position.df, f=position.df$cluster) %>% 
    lapply(reshape2::dcast,  PWM ~ y , value.var = "Freq", fill=0)  %>%
    lapply(column_to_rownames, "PWM")

p <- lapply(names(position.matrix.list), function(x, y) {
    pheatmap(y[[x]], 
             cluster_rows = TRUE, 
             cluster_cols = FALSE, 
             show_rownames = TRUE, 
             show_colnames = TRUE,
             border_color = NA,
             breaks = seq(0, 0.1, length.out =101),
             fontsize = 8, 
             fontsize_row = 8,
             fontsize_col = 8,
             main = paste("Freq for", x),
             width = 8, height = 4)
}, position.matrix.list)

ptable <- lapply(p, getElement, "gtable")

ptable2 <- list()
for(i in 1:length(ptable))
{
    ptable2[[i]] <- ptable[[i]]
    labs <- as.character(seq(-150,200, 10)) %>% lapply(function(x) c(x, rep("", 9))) %>% unlist()  %>% .[1:397]
    ptable2[[i]]$grobs[[4]]$label <- labs
}

g <- cowplot::plot_grid(plotlist = ptable2, ncol = 1,   scale = 0.95) 
ggsave(filename = file.path(out_dir, "Freq_position_PWM_per_cluster.pdf"), g, width = 16, height = 7, limitsize=FALSE)

####### 3' end shift between pA+ and pA- reads of the same PROMPT
## median_diff: median 3' end of pA+ reads - median 3' end of pA- reads, oriented on the
## PROMPT strand (positive: pA+ reads end further downstream). p-value: Wilcoxon test.
compute_TTS_shift <- function(df, s) {
    pA.positions <- rep(df$TTS, df$pA)
    nopA.positions <- rep(df$TTS, df$nopA)
    median_pos.pA <- median(pA.positions)
    median_pos.nopA <- median(nopA.positions)
    median_diff <- median_pos.pA - median_pos.nopA
    w <- wilcox.test(pA.positions, nopA.positions)
    pval <- w$p.value
    median_diff <- ifelse(s == "+", median_diff, - median_diff)
    return(c(pval, median_diff))
}

#### ZCCHC8
dhs.wilcox.ZCCHC8 <-  intersect(dhs.ZCCHC8.only, prompts.nopA.ZCCHC8)

extraction.wilcox.ZCCHC8 <- mclapply(X = dhs.wilcox.ZCCHC8,
                              FUN = function(x) {  get.signal.as.df2(dhs = x, geno = "ZCCHC8", extraction.mode = FALSE)},
                              mc.cores = 70)
names(extraction.wilcox.ZCCHC8) <- dhs.wilcox.ZCCHC8

shift.wilcox.ZCCHC8 <- mclapply(extraction.wilcox.ZCCHC8, function(x) compute_TTS_shift(x[[1]], x[[2]]), mc.cores = 70)  %>%
    lapply(function(x) data.frame(p.val = x[1], diff.median = x[2])) %>%
    mapply(mutate, ., dhs = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    mutate(padj = p.adjust(p.val, method = "BH"))

##### RBM7
dhs.wilcox.RBM7 <-  intersect(dhs.RBM7.only, prompts.nopA.RBM7)

extraction.wilcox.RBM7 <- mclapply(X = dhs.wilcox.RBM7,
                                   FUN = function(x) {  get.signal.as.df2(dhs = x, geno = "RBM7", extraction.mode = FALSE)},
                                   mc.cores = 70)
names(extraction.wilcox.RBM7) <- dhs.wilcox.RBM7

shift.wilcox.RBM7 <- mclapply(extraction.wilcox.RBM7, function(x) compute_TTS_shift(x[[1]], x[[2]]), mc.cores = 70)  %>%
    lapply(function(x) data.frame(p.val = x[1], diff.median = x[2])) %>%
    mapply(mutate, ., dhs = names(.), SIMPLIFY=FALSE) %>%
    do.call(rbind.data.frame, .) %>%
    remove_rownames() %>%
    mutate(padj = p.adjust(p.val, method = "BH"))

all.shift <- bind_rows(ZCCHC8 = shift.wilcox.ZCCHC8, 
                     RBM7 = shift.wilcox.RBM7, .id = "geno") ## %>%
    ## dplyr::filter(padj < 0.05) 

g <- ggplot(all.shift, aes(x = diff.median, colour = geno)) +
    geom_density() + 
    labs(title = "Density of median difference pA+/pA- TTS for ZCCHC8 and RBM7") +
    theme_minimal() +
    ## scale_color_manual(values = c("RBM7" = "#3D405B", "ZCCHC8" =  "#F2CC8F")) +
    scale_color_manual(values = c( "RBM7"="#b1cdc3", "ZCCHC8"="#886e97")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlim(c(-500, 1000))
ggsave(file.path(out_dir, "median_diff_density_plot.pdf"), g, width = 8, height = 5)

g <- ggplot(all.shift %>% dplyr::filter(padj <0.05), aes(x = diff.median, colour = geno)) +
    geom_density() + 
    labs(title = "Density of median difference pA+/pA- TTS for ZCCHC8 and RBM7") +
    theme_minimal() +
    scale_color_manual(values = c("RBM7" = "#3D405B", "ZCCHC8" =  "#F2CC8F")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlim(c(-500, 1000))
ggsave(file.path(out_dir, "median_diff_density_plot_signif.pdf"), g, width = 8, height = 5)

g <- ggplot(all.shift, aes(y = diff.median, x = geno, fill = geno)) +
    geom_violin() +
                                        # Point for the median
    stat_summary(fun = median, geom="point", shape=2, size=2, color="black") +
                                        # Text label for the median value
    stat_summary(fun = median, geom="text", aes(label = round(after_stat(y), 2)), 
                 vjust = -1, size = 2, color = "black") +
    ## labs(title = "Density of median difference pA+/pA- TTS for ZCCHC8 and RBM7") +
    theme_minimal() +
    scale_fill_manual(values = c( "RBM7"="#b1cdc3", "ZCCHC8"="#886e97")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(-500, 1000)) +
    coord_flip()
ggsave(file.path(out_dir, "median_diff_violin_plot.pdf"), g, width = 6, height = 3)

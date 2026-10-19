## PROMPT read 3' ends in 20 nt windows downstream of the PROMPT TSS: RBM7 / ZCCHC8
## sensitivity (DESeq2 log2FC per window) and AT content per window.
##
## Usage (from 05-PROMPT/): Rscript 06_prompt_windows.r
##
## PROMPT TSS: position covered by the most PROMPT read starts (20 nt windows; median if tied).
## Reads are counted in 60 windows of 20 nt by the position of their 3' end (0-1200 nt from
## the TSS). AT content per window and cumulative from the TSS (count and fraction).
##
## Input : $PROMPT_DIR/prompt_reads/DHS_all_reads_raw.tsv (01_prompt_reads_and_stats.r)
## Output: $PROMPT_DIR/prompt_windows/
##   at_log2FC_per_window_per_DHS.tsv   log2FC, padj, AT content, median read count per window
##   at_content_per_window_per_DHS.tsv  AT content per window
##   FC/, AT/                           heatmaps of log2FC and AT content per window
##                                      (DHSs ordered by log2FC or AT content in windows 10-30)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(parallel)
library(pheatmap)
library(data.table)
library(org.Mm.eg.db)
library(tidyverse)
library(RColorBrewer)
library(Biostrings)
library(BSgenome)

source("../config/paths.R")
source("../lib/tldr_reads.R")
source("../lib/deseq.R")
source("../04-Normalization/scale_factors.R")

out_dir <- file.path(PROMPT_DIR, "prompt_windows")
dir.create(out_dir, showWarnings=FALSE, recursive =TRUE)

## genome
genome <- Biostrings::readDNAStringSet(GENOME_FA, format = "fasta", use.names = TRUE)

dhs.gene.strand <- read.table(DHS_GENE_STRAND, sep="\t", header=TRUE)

all_reads_dt <- fread(file.path(PROMPT_DIR, "prompt_reads", "DHS_all_reads_raw.tsv"))
setnames(all_reads_dt, 
         c("V2", "V3", "V5", "V6", "V10", "V11", "V13", "V14"), 
         c("start", "end", "pA_len", "strand", "n.exon", "exon.lengths", "DHS", "rep"))

dhs.select <- dhs.gene.strand %>%
    dplyr::filter(gene_strand != "ambiguous") %>%
    .$DHS %>%
    unique

all_reads_dt[, dhs.strand := paste0(DHS, ".",strand)]
all_reads_dt <- all_reads_dt[DHS %in% dhs.select, ]

DHS.strand <- paste(dhs.gene.strand$DHS, dhs.gene.strand$gene_strand, sep=".")
all_reads_dt <- all_reads_dt[! dhs.strand %in% DHS.strand,]

find_tss_vectorized <- function(all_reads_dt) {
    ## 1. Prepare the data
    ## Ensure input is a data.table
    dt <- as.data.table(all_reads_dt %>% dplyr::select(DHS, start, end, strand, DHS))
    ## 2. Define the 20bp window anchor points based on DHS direction
    ## If DHS is '+', window starts at read start.
    ## If DHS is '-', window starts at read end - 19.
    dt[, win_start := ifelse(strand == "+", start, end - 19)]
    ## 3. Vectorized expansion
    ## Instead of looping, we expand every read into the 20 integers it covers.
    ## We use rep() for speed, creating a long-format table of every covered base.
    ## Create a vector of all start positions repeated 20 times each
    starts_rep <- rep(dt$win_start, each = 20)
    ## Create the offsets (0 to 19) repeated for every read
    offsets <- rep(0:19, times = nrow(dt))
    ## Build the expanded table
    expanded_dt <- data.table(
        DHS = rep(dt$DHS, each = 20),
        pos = starts_rep + offsets,
        strand = rep(dt$strand, each = 20)
    )
    ## 4. Calculate Coverage (Frequency Count)
    ## Count how many times each position 'pos' appears per DHS
    coverage_dt <- expanded_dt[, .N, by = .(DHS, pos, strand)]
    ## 5. Find the Peak(s) and Calculate TSS
    tss_results <- coverage_dt[, {
        max_cov <- max(N)                  ## Find max coverage value
        peak_positions <- pos[N == max_cov] ## Get all genomic positions with that max value
        .(TSS = round(median(peak_positions))) ## Calculate median of those positions
    }, by = .(DHS, strand)]
    ## 6. Merge back with metadata (to include DHSs with 0 reads)
    ## all_dhs_metadata: one row per DHS
    return(tss_results)
}

calculate_window_counts <- function(all_reads_dt, tss_results) {
    ## 1. Prepare Data
    ## We merge the reads with the calculated TSS locations
    ## Note: We assume 'all_reads_dt' has columns: DHS, start, end, strand
    dt <- merge(as.data.table(all_reads_dt %>% dplyr::select(DHS, start, end, strand, DHS, rep)), tss_results, by = c("DHS", "strand"))
    ## 2. Identify the 3' end of the reads
    ## If strand is "+", the 3' end is the 'end' coordinate.
    ## If strand is "-", the 3' end is the 'start' coordinate.
    dt[, read_3p := ifelse(strand == "+", end, start)]
    ## 3. Calculate distance from TSS
    ## We use the read's own strand to determine the direction of "downstream"
    ## (+) Strand: Distance = Read_3p - TSS
    ## (-) Strand: Distance = TSS - Read_3p
    dt[, dist_from_tss := ifelse(strand == "+", 
                                 read_3p - TSS, 
                                 TSS - read_3p)]
    ## 4. Filter and Bin
    ## We keep only reads that fall into the first 1200bp (Windows 1-60)
    ## dist_from_tss < 0 means the read 3' end is "upstream" of the TSS (discarded)
    dt_filtered <- dt[dist_from_tss >= 0 & dist_from_tss < 1200]
    ## 5. Assign Window ID (Euclidean Division)
    ## 0-19bp   -> Window 1
    ## 20-39bp -> Window 2
    dt_filtered[, window_id := (dist_from_tss %/% 20) + 1]
    ## 6. Aggregate Counts
    ## Count reads per DHS per Window
    window_counts <- dt_filtered[, .(read_count = .N), by = .(DHS, window_id, rep, TSS, strand)]
    return(window_counts)
}

calculate_window_at <- function(tss_results, genome) {
    ## 1. Prepare Data
    ## We assume tss_results has columns: DHS, TSS, strand
    dt <- as.data.table(tss_results)# %>% dplyr::filter(DHS =="chrY_1244426-1246476"))
    ## 2. Parse Chromosome from DHS ID
    ## Assumes format "chr1_123-456" or similar. 
    ## We split by "_" and take the first part.
    dt[, seqnames := tstrsplit(DHS, "_", keep = 1)]
    ## 3. Expand: Create 60 rows for every DHS
    n_dhs <- nrow(dt)
    n_wins <- 60
    win_size <- 20
    ## Expand table (repeat rows)
    expanded_dt <- dt[rep(1:n_dhs, each = n_wins)]
    ## Add window IDs
    expanded_dt[, window_id := rep(1:n_wins, times = n_dhs)]
    ## 4. Calculate Coordinates (Strand-aware)
    ## Pre-calculate offsets
    expanded_dt[, offset_start := (window_id - 1) * win_size]
    expanded_dt[, offset_end   := (window_id * win_size) - 1]
    ## (+) Strand Logic: TSS + offset,
    expanded_dt[strand == "+", `:=`(
                                   win_start = TSS + offset_start,
                                   win_end   = TSS + offset_end,
                                   window_all_start =  TSS,
                                   window_all_end = TSS + offset_end
                               )]
    ## (-) Strand Logic: TSS - offset
    ## Note: logic swaps because we move 'backward' on the genome coordinates
    expanded_dt[strand == "-", `:=`(
                                   win_end   = TSS - offset_start,
                                   win_start = TSS - offset_end,
                                   window_all_start =  TSS - offset_end,
                                   window_all_end = TSS
                               )]
    ## 5. Create GRanges
    windows_gr <- GRanges(
        seqnames = expanded_dt$seqnames,
        ranges   = IRanges(start = expanded_dt$win_start, end = expanded_dt$win_end),
        strand   = expanded_dt$strand
    )
    cumulative_at_gr <- GRanges(
        seqnames = expanded_dt$seqnames,
        ranges   = IRanges(start = expanded_dt$window_all_start, end = expanded_dt$window_all_end),
        strand   = expanded_dt$strand
    )
    ## 6. Fetch Sequences & Calculate AT
    ## This is the heavy lifting step - done in one vectorized C call
    sequences <- getSeq(genome, windows_gr)
    cumulative_sequences <- getSeq(genome, cumulative_at_gr)
    at_vals <- letterFrequency(sequences, letters = "AT", as.prob = TRUE)
    cumulative_at_vals <- letterFrequency(cumulative_sequences, letters = "AT", as.prob = TRUE)
    cumulative_at_vals.counts <- letterFrequency(cumulative_sequences, letters = "AT", as.prob = FALSE)
    widths.cumulative <- width(cumulative_sequences)
    ## 7. Format Final Output
    expanded_dt[, AT_content := as.numeric(at_vals)]
    expanded_dt[, cumulative_AT := as.numeric(cumulative_at_vals)]
    expanded_dt[, cumulative_AT.counts := as.numeric(cumulative_at_vals.counts)]
    expanded_dt[, width.cumulative := as.numeric(widths.cumulative)]
    expanded_dt[, diff_AT_20nt:= AT_content - data.table::shift(AT_content, type = "lag"), by = DHS]
    expanded_dt[, AT_40nt := frollmean(AT_content, n = 2, align = "right"), by = DHS]
    expanded_dt[, diff_AT_40nt := AT_40nt - data.table::shift(AT_40nt, n = 2, type = "lag"), by = DHS]
    expanded_dt[, AT_60nt := frollmean(AT_content, n = 3, align = "right"), by = DHS]
    expanded_dt[, diff_AT_60nt := AT_60nt - data.table::shift(AT_60nt, n = 3, type = "lag"), by = DHS]
    return(expanded_dt[, .(DHS, window_id, AT_content, cumulative_AT, diff_AT_20nt,  diff_AT_40nt, diff_AT_60nt, cumulative_AT.counts, width.cumulative)])
}

tss_results <- find_tss_vectorized(all_reads_dt)

window_at_df <- calculate_window_at(tss_results, genome)

window_counts <- calculate_window_counts(all_reads_dt, tss_results) %>%
    dplyr::filter(!grepl("ZC3H4|ZC3H18|ARS2", rep)) 

window_counts.complete <- complete(window_counts,
                                   DHS,
                                   window_id,
                                   rep,
                                   fill = list(read_count = 0, TSS = NA, strand = NA)
                                   )

count.matrix <- window_counts.complete %>%
    dplyr::select(DHS, window_id, read_count, rep) %>%
    mutate(DHS.w = paste0(DHS, ".w", window_id)) %>%
    reshape2::dcast( DHS.w ~ rep, value.var = "read_count", fill = 0) %>%
    column_to_rownames("DHS.w") 

df.meta <- colnames(count.matrix) %>%
    data.frame(sample = .) %>%
    mutate(geno = sample %>% str_replace("mES_","") %>% str_replace("_rep.",""))  %>%
    column_to_rownames("sample") %>%
    mutate(geno = factor(geno, c("WT",  "RBM7",  "ZCCHC8"))) %>%
    drop_na()
    

deseq.results <- run_deseq(count.matrix, df.meta, candidate_geno = c(  "RBM7",  "ZCCHC8"))

at.log2FC.df <- deseq.results$df_lfc_shrunk %>%
    dplyr::select(DHS, log2FoldChange , metric, padj) %>%
    mutate(geno = metric %>% str_replace("lfc_", "")) %>%
    mutate(geno = factor(geno, levels = c("ARS2", "ZC3H18",  "ZCCHC8", "RBM7")))  %>%
    separate(DHS, into = c("DHS", "window_id"), sep = "\\.w") %>%
    mutate(window_id = as.integer(window_id)) %>%
    left_join(window_at_df, by = c("DHS", "window_id")) %>%
    left_join(tss_results, by =c("DHS")) %>%
    mutate(window_id = factor(window_id, levels = 1:60))

   
##### let's try to make a heatmap of the log2FC vs window
window_counts_summary <-  window_counts %>%
    separate(rep , into = c("garb","geno", "repn"), sep = "_") %>%
    dplyr::select(-garb, -repn) %>%
    group_by(window_id, geno, DHS) %>%
    summarize(median_count= median(read_count))

fwrite(at.log2FC.df %>%
       mutate(window_id = as.integer(as.character(window_id))) %>%
       left_join(window_counts_summary, by =c("DHS", "geno", "window_id")) %>%
       mutate(median_count = ifelse(is.na(median_count), 0, median_count))
     , file.path(out_dir, "at_log2FC_per_window_per_DHS.tsv"), sep = "\t")

fwrite(window_at_df, file.path(out_dir, "at_content_per_window_per_DHS.tsv"), sep = "\t")

### try with geom raster
## Color palettes
blues_pal <- brewer.pal(n = 9, name = "Blues")
## at_pal_base <- colorRampPalette(rev(brewer.pal(n = 9, name = "Blues"))[c(1, 1:9, 9,9,9)])(100) %>% rev()
at_pal_base <- brewer.pal(n = 9, name = "Blues")[c(1,1:9,9,9)]
all.tabs <- c("FC", "AT")
for (i in 1:2){
    name.tab <- all.tabs[i]
    dir.create(file.path(out_dir, name.tab), showWarnings = FALSE)
    dir.heatmaps <- file.path(out_dir, name.tab)
    for(geno.elt in  c("RBM7", "ZCCHC8")){
        FC.window.matrix <- at.log2FC.df %>%
            dplyr::select(DHS, window_id, log2FoldChange, geno, padj) %>%
            dplyr::filter(geno == geno.elt) %>%
            mutate(logFC = log2FoldChange) %>%
            mutate(logFC = ifelse(logFC < 0, 0, logFC)) %>%
            reshape2::dcast(DHS ~ window_id, value.var = "logFC", fill = NA) %>%
            column_to_rownames("DHS") %>%
            .[rowSums(is.na(.)) < (ncol(.) - 40), ]
        dhs.filtered <- rownames(FC.window.matrix)
        window_counts.matrix <- window_counts_summary %>%
            dplyr::select(DHS, window_id, median_count, geno) %>%
            dplyr::filter(geno == geno.elt) %>%
            mutate(median_ = ifelse(median_count > 5, 5, median_count)) %>%
            reshape2::dcast(DHS ~ window_id, value.var = "median_", fill = 0) %>%
            column_to_rownames("DHS")
        AT.window.matrix <- window_at_df %>%
            dplyr::select(DHS, window_id, AT_content) %>%
            dplyr::filter(DHS %in% dhs.filtered) %>%
            mutate(AT_content = ifelse(AT_content > 0.7, 0.7, AT_content)) %>%
            mutate(AT_content = ifelse(AT_content < 0.45, 0.45, AT_content)) %>%
            reshape2::dcast(DHS ~ window_id, value.var = "AT_content") %>%
            column_to_rownames("DHS")
        list.tab <- list(FC = FC.window.matrix,
                         AT = AT.window.matrix,
                         counts = window_counts.matrix)
        tab <- list.tab[[name.tab]]
        FC.window.matrix[FC.window.matrix > 3] <- 3
        ## Order based on sum of columns 10 to 30
        order_idx <- tab[, 10:30] %>% rowSums(na.rm = TRUE) %>% order() %>% rev
        DHS.order <- rownames(tab)[order_idx]
        ## Generate custom axis labels 
        relabels <- paste0(seq(0, by = 20, length.out = 60), "-", seq(20, by = 20, length.out = 60))
        relabels[which(0:length(relabels) %% 5 != 0)] <- ""
        ## Pivot FC matrix to long format for geom_raster and apply the sorting factor
        FC.long <- FC.window.matrix %>%
            rownames_to_column("DHS") %>%
            reshape2::melt(id.vars = "DHS", variable.name = "window_id", value.name = "value") %>%
            mutate(DHS = factor(DHS, levels = DHS.order))
        p_fc <- ggplot(FC.long, aes(x = window_id, y = DHS, fill = value)) +
            geom_raster() +
            scale_fill_gradientn(colors = blues_pal, na.value = "white") +
            scale_x_discrete(labels = relabels) +
            labs(x = "Windows downstream of TSS",
                 y = "DHS") +
            theme_minimal() +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.grid = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(filename = file.path(dir.heatmaps, paste0("heatmap_log2FC_", geno.elt, "_vs_WT_across_windows_downstream_of_TSS.png")),
               plot = p_fc,
               width = 8,
               height = 10)
    }
    ## Pivot AT matrix to long format for geom_raster and apply the sorting factor
    AT.long <- AT.window.matrix %>%
        rownames_to_column("DHS") %>%
        reshape2::melt(id.vars = "DHS", variable.name = "window_id", value.name = "value") %>%
        mutate(DHS = factor(DHS, levels = DHS.order))
    p_at <- ggplot(AT.long, aes(x = window_id, y = DHS, fill = value)) +
        geom_raster() +
        scale_fill_gradientn(colors = at_pal_base, na.value = "white") +
        scale_x_discrete(labels = relabels) +
        labs(x = "Windows downstream of TSS",
             y = "DHS") +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(filename = file.path(dir.heatmaps, "heatmap_AT_content_across_windows_downstream_of_TSS.png"),
           plot = p_at,
           width = 8,
           height = 10)
}

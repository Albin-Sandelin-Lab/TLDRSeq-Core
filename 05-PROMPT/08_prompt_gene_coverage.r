## Gene and PROMPT read coverage around each DHS (heatmap, one row pair per DHS).
##
## Usage (from 05-PROMPT/): Rscript 08_prompt_gene_coverage.r
##
## For each DHS, coverage (splicing-aware, normalized by the chrM scale factors and divided by
## the number of replicates) of gene-strand reads and PROMPT reads (opposite strand), per
## genotype and pA status (pA+: tail > 15 nt). PROMPT reads extending past the gene TSS are
## removed. Positions are relative to the gene TSS (most frequent read start, 20 nt windows),
## oriented on the gene strand; +/- 2500 bp. Gene coverage is plotted as positive (blue),
## PROMPT coverage as negative (red), log10(1 + coverage), capped at 4.
## DHSs are ordered by the extent of their RBM7 pA- PROMPT signal.
##
## Input : $PROMPT_DIR/prompt_reads/DHS_all_reads_raw.tsv (01_prompt_reads_and_stats.r)
## Output: $PROMPT_DIR/prompt_gene_coverage/DHS_prompt_gene_coverage_sampled_new_order.png
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

out_dir <- file.path(PROMPT_DIR, "prompt_gene_coverage")
dir.create(out_dir, showWarnings=FALSE, recursive =TRUE)

dhs.gene.strand <- read.table(DHS_GENE_STRAND, sep="\t", header=TRUE)

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
    .[!grepl("ambiguous", .)]

all_prompt_reads_dt <- all_reads_dt2[! dhs.strand %in% DHS.strand,]
all_gene_reads_dt <- all_reads_dt2[ dhs.strand %in% DHS.strand,]
rm(all_reads_dt)
gc()

all_gene_reads_dt[, type := "gene"]
all_prompt_reads_dt[,type := "prompt"]

all_reads_dt <- rbind(all_prompt_reads_dt, all_gene_reads_dt)

all_reads_dt[, geno := tstrsplit(rep, "_", fixed = TRUE, keep = 2)] #better data.table code.
all_reads_dt[, pA_status := ifelse(pA_len > 15, "pA", "no_pA")]
all_reads_dt[, geno := paste0(geno, ":",  pA_status)]

#### Gene TSS of each DHS: position covered by the most gene-strand read starts (20 nt windows)
find_tss_ultra_safe <- function(dt) {
    dt <- as.data.table(dt)
    dt[, win_start := as.numeric(fifelse(strand == "+", start, end - 19))]

    # This processes each DHS one by one (or in small batches)
    # It never creates a massive 2.2-billion row object in memory
    tss_results <- dt[, {
        # Create a local vector for THIS group only
        # (e.g., 50 reads * 20bp = 1000 rows, very safe)
        all_pos <- as.vector(outer(win_start, 0:19, `+`))
        
        # Fast frequency count using tabulate
        # We shift by min_p to keep the index positive
        min_p <- min(all_pos)
        counts <- tabulate(all_pos - min_p + 1)
        max_n <- max(counts)
        
        # Identify the peak position(s)
        peaks <- which(counts == max_n) + min_p - 1
        
        .(TSS = floor(median(peaks)))
    }, by = .(DHS, strand)]
    
    return(tss_results)
}

tss_locations <- find_tss_ultra_safe(all_gene_reads_dt)

dhs_center.df <- all_gene_reads_dt %>%
    dplyr::select(DHS) %>%
    distinct() %>%
    separate(DHS, into = c("chrom", "start", "end"), sep = "[_-]", remove =FALSE) %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    mutate(dhs_center = round((start + end) / 2)) %>%
    dplyr::select(DHS, dhs_center)

tss_locations <- tss_locations %>%
    left_join(dhs_center.df, by="DHS") 

# Apply the function to each DHS and combine the results

df_DHS_prompt_map.gene_tss <- tss_locations %>%
    dplyr::rename(center = dhs_center, gene_start = TSS, gene_strand = strand) %>%
    mutate(offset = floor(center - gene_start))  %>%
    mutate(offset = ifelse(gene_strand =="+", -offset, offset)) 

## Reads to GRanges (one range per exon), keeping the normalization weight (scale)
reads_df_to_granges <- function(df) {
    if (is.null(df)) {
        return(GRanges())
    }
    if (nrow(df) == 0) {
        return(GRanges())
    }
    ## Initialize empty GRanges objects
    gr_unspliced <- GRanges()
    gr_spliced <- GRanges()
    ## --- Process unspliced reads (fast, no change needed) ---
    unspliced <- df[df$blockCount == 1, ]
    if (nrow(unspliced) > 0) {
        gr_unspliced <- GRanges(
            seqnames = Rle(unspliced$seqnames),
            ranges = IRanges(start = unspliced$start, end = unspliced$end),
            strand = Rle(unspliced$strand),
            scale = unspliced$scale    # <-- SURGICALLY INJECTED HERE
        )
    }
    ## --- Process spliced reads with a vectorized approach ---
    spliced <- df[df$blockCount > 1, ]
    if (nrow(spliced) > 0) {
        ## 1. Split the comma-separated strings into list-columns
        spliced_long <- spliced %>%
            mutate(
                blockSizes = str_split(blockSizes, ","),
                blockStarts = str_split(blockStarts, ",")
            )
        ## 2. Unnest the list-columns to create a "long" data frame where each
        ##    row is a single exon. This safely duplicates the 'scale' column!
        spliced_long <- spliced_long %>%
            unnest(c(blockSizes, blockStarts)) %>%
            mutate(
                blockSizes = as.numeric(blockSizes),
                blockStarts = as.numeric(blockStarts)
            )
        ## 3. Calculate exon coordinates in a vectorized way
        spliced_long <- spliced_long %>%
            mutate(
                exon_start = start + blockStarts,
                exon_end = exon_start + blockSizes - 1
            )
        ## 4. Create a single GRanges object from the long data frame.
        gr_spliced <- GRanges(
            seqnames = Rle(spliced_long$seqnames),
            ranges = IRanges(start = spliced_long$exon_start, end = spliced_long$exon_end),
            strand = Rle(spliced_long$strand),
            scale = spliced_long$scale  # <-- SURGICALLY INJECTED HERE
        )
    }
    return(c(gr_unspliced, gr_spliced))
}

### load norm factor
df_scale_factor = read.csv(SCALE_FACTORS_CSV) %>%
    mutate(Rep = sample %>% str_replace("scale_","")) %>%
    dplyr::select(scale_factor, rep = Rep) %>%
    mutate(scale = 1/scale_factor)  %>%
    dplyr::select(-scale_factor) %>%
    as.data.table()

all_reads_dt[df_scale_factor, scale := i.scale, on = "rep"]

# 1. Prepare the reads
# column names expected by reads_df_to_granges()
setnames(all_reads_dt, 
         old = c("V1", "n.exon", "exon.lengths", "V12"), 
         new = c("seqnames", "blockCount", "blockSizes", "blockStarts"), 
         skip_absent = TRUE)

# 2. Join the Reads with the DHS Metadata
# This avoids searching or filtering inside a loop.
# We only keep reads that actually map to a DHS we care about.
setkey(all_reads_dt, DHS)
setkey(df_DHS_prompt_map.gene_tss, DHS)

## 3. The Big Vectorized Calculation
## We group by (DHS, geno, type) to process everything in one shot.
raw_plot_data <- all_reads_dt[df_DHS_prompt_map.gene_tss, on = "DHS", allow.cartesian = TRUE][, {
    ## 1. Scalar assignment to prevent 'length > 1' errors
    ## Since we group by DHS, these values are constant for the group
    g_strand  <- gene_strand[1]
    g_start   <- gene_start[1]
    this_type <- type[1]
    this_dhs  <- DHS[1]
    ## --- TSS Filtering Logic ---
    keep <- TRUE
    if (this_type == "prompt") {
        if (g_strand == "+") keep <- end < g_start
        else keep <- start > g_start
    }
    ## 2. Process only if reads remain after filtering
    if (any(keep)) {
        sub_reads <- .SD[keep]
        ## Parse DHS ID for coordinates (e.g., "chr1_100-200")
        ## tstrsplit is the fast data.table way to split strings
        dhs_parts <- tstrsplit(this_dhs, "_|-")
        chrom     <- dhs_parts[[1]]
        dhs_start <- as.numeric(dhs_parts[[2]])
        dhs_end   <- as.numeric(dhs_parts[[3]])
        center    <- (dhs_start + dhs_end) / 2
        ## 3. Splicing-aware coverage
        ## reads_df_to_granges() respects blockCount/Sizes/Starts
        gr <- reads_df_to_granges(sub_reads)
        rle_list <- coverage(gr, weight = "scale")
        ## Safety check for chromosome existence in coverage object
        if (chrom %in% names(rle_list)) {
            rle_obj <- rle_list[[chrom]]
            ## Extract non-zero coverage runs (very memory efficient)
            res <- data.table(
                score = runValue(rle_obj),
                p_start = start(rle_obj),
                p_end = end(rle_obj)
            )[score > 0]
            if (nrow(res) > 0) {
                ## 4. Expand only the regions with signal (Intron-aware)
                expanded <- res[, .(
                    position = seq(p_start, p_end),
                    score = score
                ), by = p_start]
                ## Apply strand-aware relative positioning
                mult <- if (g_strand == "-") -1 else 1
                expanded[, relative_pos := (position - center) * mult]
                ## Return columns for the final combined table
                .(position = expanded$position, 
                  score = expanded$score, 
                  relative_pos = expanded$relative_pos)
            }
        }
    }
}, by = .(DHS, geno, type)]

all_dhs_ids <- unique(df_DHS_prompt_map.gene_tss$DHS)
all_types <- c("gene", "prompt")
all_genos <- c("WT",  "ARS2", "ZC3H18",  "ZC3H4",  "RBM7", "ZCCHC8") %>% paste(c(rep("pA", 5), rep("no_pA", 5)), sep=":")
#######
## Create a new score column where PROMPT coverage is negative

plot_data <- raw_plot_data %>%
    mutate(score = ifelse(geno =="WT", score/6, score/3)) %>% 
    rename(dhs_id = DHS) %>% # Change 'DHS' to 'dhs_id'
    mutate(plot_score = ifelse(type == "prompt", -log10(1 + score), log10(1 + score))) %>%
    tidyr::complete(dhs_id = all_dhs_ids, 
                    type = all_types, 
                    geno = all_genos, 
                    fill = list(relative_pos = 0, plot_score = 0))

split <- strsplit(plot_data$geno , ":", fixed = TRUE)
split.unlist <- unlist(split)
indices <- seq(1, length(split.unlist), by = 2)
geno <- split.unlist[indices]
pA_status <- split.unlist[1 + indices]
plot_data$geno <- geno
plot_data$pA_status <- pA_status
## Create a unique identifier for each row in the plot (e.g., "dhs1_Gene", "dhs1_PROMPT")
plot_data$plot_row_id <- paste(plot_data$dhs_id, plot_data$type, sep = "_")
## Define the correct order of these rows
dhs_order <- unique(all_dhs_ids)
plot_row_levels <- paste(rep(dhs_order, each = 2), c("gene", "prompt"), sep = "_")
## Convert to a factor, reversing the order so the first DHS appears at the top of the plot
plot_data$plot_row_id <- factor(plot_data$plot_row_id, levels = rev(plot_row_levels))
## Convert geno to factor for facet ordering
plot_data$geno <- factor(plot_data$geno, levels = c("WT", "RBM7", "ZCCHC8"))
plot_data <- plot_data %>% drop_na
######

## same but with data.table functions

plot_data.sample <- plot_data    %>% dplyr::filter(dhs_id %in% unique(dhs_id))

# --- 4. Generate the Plot ---

max_dist <- 2500

plot_data_filtered <- plot_data.sample %>%
     left_join(df_DHS_prompt_map.gene_tss, by=c("dhs_id"="DHS"), multiple = "first") %>%
    mutate(relative_pos = relative_pos - offset) 

plot_data_filtered <- plot_data_filtered %>%
  filter(abs(relative_pos) <= max_dist)
### reorder the rows again
dhs_order <- plot_data_filtered %>%
    dplyr::filter(plot_score < 0, geno =="RBM7", pA_status == "no_pA", type =="prompt") %>%
    group_by(dhs_id) %>%
    mutate(min.score.dhs = min(plot_score)) %>%
    ungroup %>%
    dplyr::filter(plot_score < (min.score.dhs/10)) %>%
    group_by(dhs_id) %>% 
    slice_max(relative_pos) %>%
    arrange(desc(relative_pos)) %>% 
    .$dhs_id

dhs_order.end <- plot_data_filtered %>%
    dplyr::filter(dhs_id %in% setdiff(unique(plot_data_filtered$dhs_id), dhs_order)) %>%
    dplyr::filter(plot_score > 0, geno =="RBM7", pA_status == "no_pA", type =="gene") %>%
    group_by(dhs_id) %>%
    slice_max(plot_score, with_ties=FALSE ) %>%
    arrange(desc(plot_score)) %>%
    .$dhs_id

plot_row_levels <- paste(rep(c( dhs_order), each = 2), c("gene", "prompt"), sep = "_")

## Convert to a factor, reversing the order so the first DHS appears at the top of the plot
plot_data_filtered$plot_row_id2 <- factor(plot_data_filtered$plot_row_id, levels = rev(plot_row_levels)) 
plot_data_filtered <- plot_data_filtered %>%
    dplyr::filter(!is.na(plot_row_id2))

plot_data_filtered$geno <- factor(plot_data_filtered$geno, levels = c("WT",    "RBM7", "ZCCHC8"))

## Get the range of scores to set symmetrical color limits for the legend
max_abs_score <- 4 #max(abs(plot_data_filtered$plot_score), na.rm = TRUE)
plot_data_filtered$plot_score <- ifelse(abs(plot_data_filtered$plot_score ) > 4, 
                                        sign(plot_data_filtered$plot_score) * 4, 
                                        plot_data_filtered$plot_score)

rdbu <- rev(brewer.pal(n = 11, name = "RdBu"))

p <- ggplot(plot_data_filtered, aes(x = relative_pos, y = plot_row_id2)) +
    geom_tile(aes(fill = plot_score)) +
    scale_fill_gradientn(colours=c(rdbu[2], rdbu[3], rdbu[4], "white" ,"white", rdbu[8], rdbu[9], rdbu[10]),
                          values=scales::rescale(c(-4,  -1.5, -0.3, -0.05,  0.05, 0.3, 1.5, 4)),
                          breaks=c(-4,  -1.5, 0, 1.5, 4),
                          limits=c(-4,4)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
    facet_grid(rows = vars(pA_status), cols = vars(geno)) +
    labs(
        x = "Position relative to DHS center (bp)",
        y = "DHS Sites (Gene/PROMPT pairs)",
        title = "Gene and PROMPT Read Coverage at DHS Sites",
        subtitle = "Genes are positive (blue), PROMPTs are negative (red)"
    ) +
    scale_x_continuous(labels = scales::comma, breaks = c(-2000, 0, 2000)) +
    theme_minimal() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "grey90", color = "grey90"),
        strip.text = element_text(face = "bold", size = 10)
    )

ggsave(filename = file.path(out_dir, "DHS_prompt_gene_coverage_sampled_new_order.png"),
       plot = p,
       width = 10,
       height = 12)

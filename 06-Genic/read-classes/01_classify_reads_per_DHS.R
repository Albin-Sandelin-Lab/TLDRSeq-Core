## Assign TLDR-Seq reads to DHSs (by their 5' end) and classify them by their 3' end
## relative to the features of the associated gene.
##
## Usage (from 06-Genic/read-classes/): Rscript 01_classify_reads_per_DHS.R <dataset> [n_cores]
##   dataset : TLDRSeq_sup, TLDRSeq_nuclear_sup or TLDRSeq_Z1_MTR4_sup
##   n_cores : samples processed in parallel (default 8). Use n_cores = 1 if too much load
##             (memory: the reads of each sample overlapping the DHSs are loaded at once).
##
## Input : $RESULTS_DIR/<dataset>/concatenate_per_replicate/<sample>.bam   (01-TLDRSeq-alignment)
##         $RESULTS_DIR/<dataset>/separate_pA_reads/reads_pA_length.tsv     (01-TLDRSeq-alignment)
##         DHS gene features: DHS_GENE_MAP_EXONS, DHS_TSS_CURATED, DHS_TTS_CURATED (config/paths.R)
## Output: $GENIC_DIR/<dataset>/attribute_reads/df_unified_results_classified.rds
##         one row per read on the gene strand: DHS, gene features, exon_count, polyA_length,
##         label_3prime, polyA_status, label_3prime_pA, sample,
##         depth.dhs (reads assigned to DHSs in the sample, both strands)
##
## Read categories (label_3prime), from the read 3' end:
##   upstream              ends upstream of the curated TSS
##   typical_TTS           within +/- 100 bp of the curated TTS
##   first_exon            unspliced, before the TTS, ends in the first exon
##   premature_termination unspliced, before the TTS, ends in the first intron
##   gradual_degradation   before the TTS, ends after the first intron and before the CDS end
##   over_elongation       more than 100 bp after the TTS
##   Others                everything else (3' UTR, spliced reads ending in the first exon
##                         or intron)
## Only reads on the strand of the gene are kept. DHSs of single-exon genes are discarded.
## pA+: poly(A) tail > 15 nt.

library(GenomicAlignments)
library(GenomicRanges)
library(tidyverse)
library(rtracklayer)
library(data.table)

source("../../config/paths.R")

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
n_cores <- ifelse(length(args) >= 2, as.integer(args[2]), 8)

# --- 1. SETTINGS & PATHS ---
bam_dir <- file.path(RESULTS_DIR, dataset, "concatenate_per_replicate")
pa_dir  <- file.path(RESULTS_DIR, dataset, "separate_pA_reads")
out_dir <- file.path(GENIC_DIR, dataset, "attribute_reads")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load Curated Annotations
gr_tts   <- readRDS(DHS_TTS_CURATED)
gr_tss   <- readRDS(DHS_TSS_CURATED)
df_map   <- readRDS(DHS_GENE_MAP_EXONS)

# Load PolyA Data
pa_file <- file.path(pa_dir, "reads_pA_length.tsv")
pa_data <- fread(pa_file) %>%
    as_tibble() %>%
    setNames(c("ID", "Value")) %>%
    mutate(Value = str_replace(Value, "pt:i:", "") %>% as.numeric())

# Create DHS GRanges
dhs_coords <- str_match(df_map$DHS, "(chr[0-9XYM]+)_([0-9]+)-([0-9]+)")
dhs_gr <- GRanges(
  seqnames = dhs_coords[,2],
  ranges = IRanges(as.numeric(dhs_coords[,3]), as.numeric(dhs_coords[,4])),
  name = df_map$DHS
)
## Regions to load from the BAM files (a read whose 5' end is in a DHS overlaps it)
dhs_gr_reduced <- reduce(dhs_gr, ignore.strand = TRUE)


# --- 2. THE MASTER PROCESSING FUNCTION ---
process_replicate <- function(sample_id) {
    message(paste(">>> Processing Sample:", sample_id))
    bam_file <- file.path(bam_dir, paste0(sample_id, ".bam"))
    if(!file.exists(bam_file)) {
        warning(paste("BAM not found for", sample_id)); return(NULL)
    }

    # Load reads (primary alignments overlapping the DHSs only)
    param <- ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE),
                          which = dhs_gr_reduced)
    ga <- readGAlignments(bam_file, use.names = TRUE, param = param)
    ## a read overlapping several regions of 'which' is returned once per region: keep one copy
    ga <- ga[!duplicated(paste(names(ga), seqnames(ga), start(ga), strand(ga), cigar(ga)))]
    reads_grl <- granges(ga, use.mcols = TRUE)

    # A. Get junction/exon counts
    mcols(reads_grl)$exon_count <- njunc(ga) + 1

    # B. Global PolyA Match
    mcols(reads_grl)$polyA_length <- pa_data$Value[match(names(reads_grl), pa_data$ID)]
    mcols(reads_grl)$polyA_length <- ifelse(is.na(mcols(reads_grl)$polyA_length), 0, mcols(reads_grl)$polyA_length)

    # C. Global DHS Filter (5' end anchoring)
    tss_reads <- resize(reads_grl, width = 1, fix = "start")
    hits <- findOverlaps(tss_reads, dhs_gr)
    valid_reads <- reads_grl[queryHits(hits)]
    mcols(valid_reads)$DHS <- dhs_gr$name[subjectHits(hits)]
    depth.DHS <- length(valid_reads)

    # D. High-Speed Classification
    df_results <- as.data.frame(valid_reads) %>%
        as_tibble() %>%
        left_join(df_map, by = "DHS") %>%

        ## Keep reads on the gene strand (strand.x: read, strand.y: gene)
        filter(strand.x == strand.y) %>%

        ## -----------------------------------------------------------
        ## RULE 0: STRICTLY DISCARD SINGLE EXON GENES
        ## -----------------------------------------------------------
        filter(!is.na(second_exon_start)) %>%

        mutate(
            ref_tts = start(gr_tts)[match(DHS, gr_tts$DHS)],
            ref_tss = start(gr_tss)[match(DHS, gr_tss$DHS)]
        ) %>%
        mutate(
            read_3p = ifelse(strand.x == "+", end, start),
            dist_to_tts = ifelse(strand.x == "+", read_3p - ref_tts, ref_tts - read_3p)
        ) %>%
        mutate(
            ## -----------------------------------------------------------
            ## DEFINE GEOGRAPHIC BOUNDARIES & FACTS
            ## -----------------------------------------------------------
            unspliced = (exon_count == 1),

            # Distance from annotated TTS
            label_tts_fact = case_when(
                abs(dist_to_tts) <= 100 ~ "typical_TTS",
                dist_to_tts > 100       ~ "overlong_TTS",
                dist_to_tts < -100      ~ "before_TTS"
            ),

            # Position relative to gene features
            upstream_of_tss = ifelse(strand.x == "+", read_3p < ref_tss, read_3p > ref_tss),
            is_before_cds   = ifelse(strand.x == "+", read_3p < cdsEnd, read_3p > cdsStart),

            # Intron/Exon 1 Boundaries
            is_in_first_exon = ifelse(strand.x == "+",
                                      read_3p <= first_exon_end,
                                      read_3p >= first_exon_end),

            is_in_first_intron = ifelse(strand.x == "+",
                                        read_3p > first_exon_end & read_3p < second_exon_start,
                                        read_3p < first_exon_end & read_3p > second_exon_start),

            is_after_first_intron = ifelse(strand.x == "+",
                                           read_3p >= second_exon_start,
                                           read_3p <= second_exon_start)
        ) %>%
        mutate(
            ## -----------------------------------------------------------
            ## BIOLOGICAL CATEGORIZATION (The 4 Core Categories + Others)
            ## -----------------------------------------------------------
            label_3prime = case_when(

                # Intercept Upstream reads immediately (so they don't get mislabeled as 1st Exon)
                upstream_of_tss ~ "upstream",

                ## CATEGORY 1: Typical TTS (Within +/- 100bp window)
                label_tts_fact == "typical_TTS" ~ "typical_TTS",

                ## CATEGORY 2: First Exon (Strictly unspliced, ends in Exon 1)
                label_tts_fact == "before_TTS" & unspliced & is_in_first_exon ~ "first_exon",

                ## CATEGORY 3: Premature Termination (Strictly unspliced, ends in Intron 1)
                label_tts_fact == "before_TTS" & unspliced & is_in_first_intron ~ "premature_termination",

                ## CATEGORY 4: Gradual Degradation (Ends at/after Exon 2 AND before CDS stop)
                # Note: Splicing agnostic.
                label_tts_fact == "before_TTS" & is_after_first_intron & is_before_cds ~ "gradual_degradation",

                ## Over-elongation (ends more than 100 bp after the TTS)
                label_tts_fact == "overlong_TTS" ~ "over_elongation",

                ## CATCH-ALL: Others
                # Absorbs: 3' UTR degradation, and
                # the "weird" artifact reads (spliced reads that terminate in Exon 1 or Intron 1)
                TRUE ~ "Others"
            ),

            # Add PolyA flags
            polyA_status = ifelse(polyA_length > 15, "pA+", "pA-"),
            label_3prime_pA = paste0(label_3prime, "_", polyA_status),
            sample = sample_id,
            depth.dhs = depth.DHS
        )

    return(df_results)
}

# --- 3. RUN PIPELINE ---
## All per-replicate BAMs of the dataset (01-TLDRSeq-alignment/03_merge_per_sample.sh)
sample_list <- list.files(bam_dir, pattern = "\\.bam$") %>% str_replace("\\.bam$", "")

# Map over samples and combine into one master dataframe
final_classified_data <- sample_list %>%
    parallel::mclapply(process_replicate, mc.cores = n_cores) %>%
    bind_rows()

# Save final results
saveRDS(final_classified_data, file.path(out_dir, "df_unified_results_classified.rds"))

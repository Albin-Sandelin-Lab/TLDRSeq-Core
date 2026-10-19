## Shared paths for the R scripts (same values as config/paths.sh).
## Edit this file to run on another system.
## R scripts are run from their step folder and source "../config/paths.R".

DATA_DIR    <- "/maps/projects/sandelin_main/data/degradation_project"
RESULTS_DIR <- "/maps/projects/sandelin_main/people/fzh976/degradation_project/results"

## Genome and annotations (mm39, chr-prefixed names)
GENOME_FA   <- "/maps/projects/sandelin_main/data/shared_data/genomes/mm39/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GENCODE_GTF <- "/maps/projects/sandelin_main/data/shared_data/genomes/mm39/gencode.vM35.primary_assembly.annotation.gtf"

#### Inputs produced by the DHS annotation / read assignment (03-Annotation, not in this repository yet)
ANNOTATION_DIR <- "/maps/projects/sandelin_main/people/fzh976/degradation_project/scripts/dhs_charac/depth_around_dhs/Rscripts"

## DHS annotation: overlaps with RefSeq select / GENCODE genes, TSS, TTS, MGI ncRNAs
DHS_ANNOTATION_RDS <- file.path(ANNOTATION_DIR, "20251028_annotation_of_DHS_with_everything.rds")
## Strand of the gene associated with each DHS (gene_strand, gencode.intersect)
DHS_GENE_STRAND    <- file.path(RESULTS_DIR, "dhs_charac/annotate_DHS/DHS_gene_strand.tsv")
## RefSeq select annotation (GTF)
REFSEQ_SELECT_GTF  <- file.path(ANNOTATION_DIR, "mouse.mm39.refseq.select.gtf")
## TLDR-Seq reads starting at each DHS, one folder per DHS and replicate
DIR_READS_PER_DHS  <- file.path(RESULTS_DIR, "dhs_charac/reads_start_from_each_DHS_20250331_run_after_DHS_refine")
## chrM read depth per replicate, used for normalization (scale factors)
CHRM_DEPTH         <- file.path(RESULTS_DIR, "dhs_charac/get_bed12_per_reads/bed12_primary_alignment_chrM_per_replicate/read_depth_chrM.txt")
## chrM scale factor per replicate (04-Normalization/01_write_scale_factors.R)
SCALE_FACTORS_CSV  <- file.path(RESULTS_DIR, "04-Normalization", "scale_factor_per_rep.csv")

## Gene features of DHSs uniquely associated with one gene (2025-11-18): curated TSS and TTS
## (peak calling on TLDR-Seq read ends), first exon / second exon, CDS
DHS_GENE_MAP_EXONS <- file.path(ANNOTATION_DIR, "20251118_df_DHS_gene_unique_map_with_exon_info.rds")
DHS_TSS_CURATED    <- file.path(ANNOTATION_DIR, "20251118_grange_tss_curated_by_peakcalling.rds")
DHS_TTS_CURATED    <- file.path(ANNOTATION_DIR, "20251118_grange_tts_curated_by_peakcalling.rds")

#### Enhancer DHSs (07-Enhancer)
## Enhancer DHS annotation (TSS of the + and - strand transcripts, DHS width) and public bigWig
## tracks lifted to mm39 (ATAC, H3K27ac, H3K27me2/3, H3K4me1, STARR-seq, DNase; GEO)
ENHANCER_DATA_DIR       <- "/maps/projects/sandelin_main/data/histones_mm39"
ENHANCER_ANNOTATION_RDS <- file.path(ENHANCER_DATA_DIR, "20260323_df_enhancer_DHS_with_properties.rds")
## H3K27ac and H3K4me3 tracks (mm39)
HISTONE_MARKS_DIR       <- file.path(DATA_DIR, "dhs_charac/histone_marks")
## ReMap mESC ChIP-seq, WT samples pooled per TF (mm39), and the area under each track
REMAP_DIR               <- file.path(DATA_DIR, "bw_remap")
REMAP_BW_DIR            <- file.path(REMAP_DIR, "pooled_WT_bigWig")
REMAP_NORM_FACTORS      <- file.path(REMAP_DIR, "TF_normalization_factors.txt")
## GSE89209 ChIP-seq (mESC and AntNPC: H3K27ac, p300, H3K27me3), bigWigs of raw read coverage (mm39)
CHIP_GSE89209_DIR       <- file.path(RESULTS_DIR, "diverse/ChIP_GSE89209")

## Output of each step
PROMPT_DIR   <- file.path(RESULTS_DIR, "05-PROMPT")
GENIC_DIR    <- file.path(RESULTS_DIR, "06-Genic")
ENHANCER_DIR <- file.path(RESULTS_DIR, "07-Enhancer")
QC_DIR       <- file.path(RESULTS_DIR, "08-QC")

## Genic read classes of TLDRSeq_sup (06-Genic): reads per DHS, replicate and 3' end label
## combination (06-Genic/read-classes/label_reads_TLDRSeq_sup.R), and frequency of each
## category per DHS and genotype, from scaled counts pooled over replicates
## (06-Genic/proportions/02_build_freq_table.R)
TTS_PROP_CATEGORIES_RDS <- file.path(GENIC_DIR, "TLDRSeq_sup", "read_labels", "TTS_prop_categories.rds")
GENIC_FREQ_RDS          <- file.path(GENIC_DIR, "TLDRSeq_sup", "proportions", "genic_freq_from_pooled_scaled_counts.rds")

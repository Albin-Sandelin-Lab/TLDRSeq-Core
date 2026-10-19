# 07 — Enhancers: external data

Download and processing of the public datasets used in 07-Enhancer. Run once, before the scripts in `Enhancer-analysis/`.
Paths are set in [`../../../config/paths.sh`](../../../config/paths.sh) and [`../../../config/paths.R`](../../../config/paths.R).

## Steps

| Script | Does | Output |
|---|---|---|
| `01_remap_download.r` | ReMap 2022 mouse ESC bigWig archive (mm10), Recherche Data Gouv doi:10.57745/GEIHVN | `$REMAP_DIR/*.tar.part.a?` |
| `02_remap_extract_lift.sh` | extract the archive; lift every bigWig over to mm39 (CrossMap) | `$REMAP_DIR/bigwig_files_mm39/` |
| `03_remap_select_pool.r` | keep experiments in control / wild-type conditions; sum the coverage of the experiments of each TF | `$REMAP_DIR/cleaned_mapping_table.txt`, `$REMAP_DIR/pooled_WT_bigWig/` |
| `04_track_areas.r` | area under each track (ReMap TFs and the tracks below), for normalization | `$REMAP_DIR/TF_normalization_factors.txt` |
| `05_chip_GSE89209_download.sh` | GSE89209 ChIP-seq (mESC, AntNPC: H3K27ac, p300, H3K27me3, input) from SRA | `$CHIP_GSE89209_FASTQ/` |
| `06_chip_GSE89209_map.sh` | bowtie2 (mm39), mapped reads, bigWig of the read coverage (10 bp bins, no normalization) | `$CHIP_GSE89209_DIR/` |
| `07_nascent_download_lift.sh <GSE>` | nascent RNA bigWigs from GEO (GSE142312 GRO-seq, GSE115774 4sU-seq), lifted from mm9 to mm39 | `$DATA_DIR/<GSE>_nascent/` |

Run order: 01 → 02 → 03 → 04 (04 also needs the tracks below); 05 → 06; 07 for each GSE.

## Other tracks (downloaded from GEO, mm39; accessions in the paper)

| File | Assay |
|---|---|
| `ENHANCER_DATA_DIR/ATAC_merged.bw` | ATAC-seq |
| `ENHANCER_DATA_DIR/H3K27ac_merged.bw` | H3K27ac ChIP-seq |
| `ENHANCER_DATA_DIR/H3K27me2_merged.bw` | H3K27me2 ChIP-seq |
| `ENHANCER_DATA_DIR/H3K27me3_merged.bw` | H3K27me3 ChIP-seq |
| `ENHANCER_DATA_DIR/H3K4me1_merged.bw` | H3K4me1 ChIP-seq |
| `ENHANCER_DATA_DIR/STARRseq_merged.bw` | STARR-seq |
| `ENHANCER_DATA_DIR/Dnase1.bw` | DNase-seq |
| `HISTONE_MARKS_DIR/H3K27ac-mm_norm_mm39.bw` | H3K27ac ChIP-seq |
| `HISTONE_MARKS_DIR/H3K4me3-mm_norm_mm39.bw` | H3K4me3 ChIP-seq |

## Requirements

R: dataverse, tidyverse, GenomicRanges, rtracklayer.
Command line: CrossMap (conda environment `$CROSSMAP_ENV`), sratoolkit, Entrez Direct, bowtie2,
samtools, deepTools, wget.

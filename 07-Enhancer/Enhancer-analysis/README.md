# 07 — Enhancers: analysis

Enhancer DHSs with a TSS on both strands. TSS center = middle of the + and − strand TSSs;
TSS distance = TSS+ − TSS− (negative: the two transcripts overlap). DHSs are sorted by TSS
distance and split into 9 quantiles.

## Run

```bash
./run_all.sh
```

Paths are set in [`../../config/paths.R`](../../config/paths.R); results go to `$RESULTS_DIR/07-Enhancer/`.

## Inputs

Public datasets are downloaded and processed by [`data_preparation/`](data_preparation/README.md) (run first).

| Input | From |
|---|---|
| enhancer DHS annotation (`ENHANCER_ANNOTATION_RDS`) | 03-Annotation |
| TLDR-Seq reads per DHS and replicate (`DIR_READS_PER_DHS`), chrM scale factors (`SCALE_FACTORS_CSV`) | 03-Annotation, 04-Normalization |
| bigWig tracks (mm39): ATAC, H3K27ac, H3K27me2/3, H3K4me1, H3K4me3, STARR-seq, DNase (`ENHANCER_DATA_DIR`, `HISTONE_MARKS_DIR`) | GEO |
| GSE89209 ChIP-seq, mESC and AntNPC (H3K27ac, p300): bigWigs of raw read coverage (`CHIP_GSE89209_DIR`) | GEO |
| nascent RNA bigWigs lifted to mm39: GSE142312 GRO-seq, GSE115774 4sU-seq (`DATA_DIR`) | GEO |
| ReMap mESC ChIP-seq, WT samples pooled per TF, lifted to mm39 (`REMAP_BW_DIR`); area under each track (`REMAP_NORM_FACTORS`) | ReMap |

## Steps

| Script | Does | Output (in `$RESULTS_DIR/07-Enhancer/`) |
|---|---|---|
| `01_enhancer_coverage_and_chip.r` | RBM7 TLDR-Seq read coverage (whole reads per strand; read 5′ ends, both strands) around the TSS center; mean signal of each bigWig track per position and TSS distance quantile, anchored on the TSS center and on each TSS; median STARR-seq signal around the TSS center of convergent DHSs; RBM7 mean read length and TSS distance per quantile | `coverage_and_chip/` |
| `02_enhancer_chip_logFC.r` | H3K27ac and p300 (GSE89209) around the TSS center, normalized by the area under each track; log2FC AntNPC vs mESC per DHS; fraction of DHSs with log2FC > 0 per TSS distance quantile | `chip_logFC/` |
| `03_nascent_tss_edges.r <dataset>` | intergenic enhancers: TSS of the + and − strand transcripts from nascent RNA (GSE142312 GRO-seq, ESC vs NPC; GSE115774 4sU-seq, Day0 vs Day3; replicates pooled); fraction of enhancers switching from convergent to divergent (coverage threshold set per dataset: 100 GRO-seq, 200 4sU-seq) | `nascent_tss_edges/<dataset>/` |
| `04_enhancer_reads_sensitivity.r` | enhancer TLDR-Seq reads: DESeq2 RBM7 / ZCCHC8 vs WT per DHS and strand; exonic length and exon number of upregulated DHSs (pA+ / pA−); pA+ read fraction per DHS and per k-means log2FC interval (ANOVA); log2FC per read length category | `reads_sensitivity/` |
| `05_enhancer_tss_logo.r` | sequence logos from −100 to +20 nt around the + and − strand TSSs of the enhancers, per TSS distance quantile | `tss_logo/` |

The heat map of the bigWig signal per TSS distance quantile (mean signal across the enhancers of a
quantile, scaled to min = 0, max = 1 per track) is made from `metadata_summary_dt_tss.tsv` and
`tss_center_summary.tsv` outside this pipeline.

## Requirements

R: tidyverse, data.table, GenomicRanges, rtracklayer, RColorBrewer, scales, ggridges, patchwork,
pbapply, plyranges, DESeq2, apeglm, classInt (needs e1071), Biostrings, BSgenome, ggseqlogo, ggpubr.

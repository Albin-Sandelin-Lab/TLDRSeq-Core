# 06 — Genic: read classes

TLDR-Seq reads starting at a DHS uniquely associated with one gene, on the gene strand,
classified by the position of their 3′ end relative to the gene features.

## Run

```bash
./run_all.sh [n_cores]
```

Paths are set in [`../../config/paths.R`](../../config/paths.R); results go to `$RESULTS_DIR/06-Genic/<dataset>/`.

## Inputs from other steps

| Input | From |
|---|---|
| per-replicate BAM files, poly(A) tail lengths | 01-TLDRSeq-alignment |
| TLDR-Seq reads per DHS and replicate (`DIR_READS_PER_DHS`) | 03-Annotation |
| DHS gene features: first / second exon, CDS, curated TSS and TTS (`DHS_GENE_MAP_EXONS`, `DHS_TSS_CURATED`, `DHS_TTS_CURATED`) | 03-Annotation, 06-Genic/peak-calling |
| chrM read depth per replicate (`CHRM_DEPTH`, `SCALE_FACTORS_CSV`) | 04-Normalization |

## Read categories

From the read 3′ end (pA+ = poly(A) tail > 15 nt):

| Category | Read 3′ end |
|---|---|
| `upstream` | upstream of the curated TSS |
| `typical_TTS` | at the curated TTS |
| `first_exon` | unspliced read, before the TTS, in the first exon |
| `premature_termination` | unspliced read, before the TTS, in the first intron |
| `gradual_degradation` | before the TTS, after the first intron and before the CDS end |
| `over_elongation` | after the TTS |

The other reads are grouped together.

## Steps

| Step | Script | Dataset | Does | Output (in `$RESULTS_DIR/06-Genic/`) |
|---|---|---|---|---|
| 1 | `label_reads_TLDRSeq_sup.R` | `TLDRSeq_sup` | label the 3′ end of the reads of each DHS (`label_reads_by_3prime.R`); reads per DHS, replicate and label combination | `TLDRSeq_sup/read_labels/` |
| 2 | [`../proportions/02_build_freq_table.R`](../proportions) | `TLDRSeq_sup` | category of each read; chrM-scaled counts summed over replicates; frequency of each category per DHS and genotype | `TLDRSeq_sup/proportions/` |
| 3 | `01_classify_reads_per_DHS.R <dataset>` | `TLDRSeq_nuclear_sup`, `TLDRSeq_Z1_MTR4_sup` | assign reads to DHSs by their 5′ end and classify them | `<dataset>/attribute_reads/` |
| 4 | `02_library_size_factors.r <dataset>` | `TLDRSeq_nuclear_sup`, `TLDRSeq_Z1_MTR4_sup` | bambu gene counts: library depth (for CPM) and chrM scale factors (for DESeq2) | `<dataset>/library_size_factors/` |
| 5 | `03_read_categories_heatmap.r` | `TLDRSeq_sup` | heatmap of the WT category frequencies; fraction of over-elongated reads in WT | `TLDRSeq_sup/read_categories/` |
| 6 | `04_cyto_nuc_ranks.r` | `TLDRSeq_nuclear_sup` | DHS × category ranked by cytoplasmic / nuclear ratio of mean CPM (no scaling between fractions: no spike-in, no mitochondrial RNA in the nuclear fraction); rank heatmaps and rank tables | `TLDRSeq_nuclear_sup/cyto_nuc_ranks/` |
| 7 | `05_read_categories_sensitivity.r` | `TLDRSeq_sup`, `TLDRSeq_Z1_MTR4_sup` | pA+ fraction per category (WT, RBM7, ZCCHC8); cytoplasmic / nuclear rank per category; DESeq2 log2FC per category for RBM7 / ZCCHC8 and MTR4 / Z1-FL; log2FC of first exonic pA− reads of cytoplasmic vs nuclear DHSs | `read_categories_sensitivity/` |
| 8 | `06_AT_model_polynomial.r` | `TLDRSeq_sup` | unspliced first exonic (nuclear DHSs only) and first intronic gene reads counted in 20 nt windows from the gene TSS; DESeq2 log2FC per window (RBM7, ZCCHC8); with the PROMPT windows (05-PROMPT), one model `log2FC ~ poly(log10(cumulative AT count), 2) + poly(log10(length), 2)`; length-corrected log2FC vs AT count, per genotype | `TLDRSeq_sup/AT_model_polynomial/` |

`06_AT_model_polynomial.r` also uses outputs of 05-PROMPT (`prompt_reads/DHS_all_reads_raw.tsv`,
`prompt_windows/at_log2FC_per_window_per_DHS.tsv`): run 05-PROMPT first.

## Requirements

R: tidyverse, data.table, GenomicRanges, GenomicAlignments, GenomicFeatures, plyranges,
rtracklayer, Biostrings, BSgenome, bambu, DESeq2, apeglm, pheatmap, RColorBrewer, reshape2, ggExtra,
pbapply.

# 05 — PROMPTs

Promoter upstream transcripts (PROMPTs): TLDR-Seq reads starting at a gene-associated
DHS on the strand opposite to the gene. Analyses use the `TLDRSeq_sup` dataset
(WT, RBM7 and ZCCHC8)

## Run

```bash
./run_all.sh
```

Paths are set in [`../config/paths.R`](../config/paths.R); results go to `$RESULTS_DIR/05-PROMPT/`.

## Inputs from other steps

| Input | From |
|---|---|
| DHS annotation (`DHS_ANNOTATION_RDS`, `DHS_GENE_STRAND`, `REFSEQ_SELECT_GTF`) | 03-Annotation |
| TLDR-Seq reads per DHS and replicate (`DIR_READS_PER_DHS`) | 03-Annotation |
| chrM scale factors (`CHRM_DEPTH`, `SCALE_FACTORS_CSV`) | 04-Normalization |

## Steps

| Script | Does | Output (in `$RESULTS_DIR/05-PROMPT/`) |
|---|---|---|
| `01_prompt_reads_and_stats.r` | select DHSs uniquely associated with one gene TSS; collect their reads; keep PROMPT reads; read counts (normalized by chrM scale factors) and read statistics (length, exons, pA length) per DHS, genotype and pA status (pA+ = tail > 15 nt) | `prompt_reads/` |
| `02_prompt_sensitivity.r` | DESeq2 RBM7 / ZCCHC8 vs WT; upregulated PROMPTs (log2FC > 1, padj < 0.05); pA fraction, read length, exon number; sensitivity heatmaps | `sensitivity/` |
| `03_prompt_vs_pmt_length.r` | PROMPT vs first exon vs first intron (PMT) vs gene reads, by read length category: DESeq2, pA fraction, upregulated fraction | `prompt_vs_pmt_length/` |
| `04_prompt_length_categories.r` | DESeq2 WT / RBM7 / ZCCHC8 per read length category; fraction of upregulated PROMPTs (pA+ / pA−) and PROMPT log2FC per length category | `prompt_length_categories/` |
| `05_prompt_pA_motifs.r` | pA+ / pA− PROMPTs in WT, RBM7, ZCCHC8: 3′ end peaks, sequence logos, XSTREME motifs and positions of motif hits, 3′ end shift between pA+ and pA− reads. PROMPT reads overlapping the gene body are excluded here | `prompt_pA_motifs/` |
| `06_prompt_windows.r` | PROMPT TSS; PROMPT reads counted in 20 nt windows by their 3′ end (0–1200 nt); DESeq2 log2FC (RBM7, ZCCHC8) and AT content per window; heatmaps | `prompt_windows/` |
| `07_prompt_AT_model.r` | RBM7 and ZCCHC8 log2FC per window vs cumulative AT count from the PROMPT TSS, corrected for length (linear model) | `AT_model/` |
| `08_prompt_gene_coverage.r` | heatmap of gene (positive) and PROMPT (negative) read coverage around each DHS, per genotype and pA status | `prompt_gene_coverage/` |

## Requirements

R: tidyverse, data.table, GenomicRanges, GenomicFeatures, rtracklayer, plyranges, Biostrings,
BSgenome, DESeq2, apeglm, pbapply, pheatmap, RColorBrewer, classInt, ggVennDiagram, ggExtra,
ggpubr, ggseqlogo, cowplot, ggpmisc, hexbin, org.Mm.eg.db.
MEME suite (`xstreme`).

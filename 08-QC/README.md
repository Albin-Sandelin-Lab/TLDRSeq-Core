# 08 — Quality control

Library size and PCA of the TLDR-Seq read counts per DHS, for `TLDRSeq_sup` (WT, RBM7, ZCCHC8).

## Run

```bash
./run_all.sh [n_cores]
```

Paths are set in [`../config/paths.R`](../config/paths.R); results go to `$RESULTS_DIR/08-QC/`.

## Inputs from other steps

| Input | From |
|---|---|
| per-replicate BAM files and their index | 01-TLDRSeq-alignment |
| DHS annotation (`DHS_ANNOTATION_RDS`), TLDR-Seq reads per DHS and replicate (`DIR_READS_PER_DHS`) | 03-Annotation |

## Steps

| Script | Does | Output (in `$RESULTS_DIR/08-QC/`) |
|---|---|---|
| `01_library_size.r` | number of mapped reads of each BAM (from its index); plot for WT, RBM7, ZCCHC8 | `library_size/` |
| `02_DHS_read_counts_PCA.r [n_cores]` | reads per replicate, DHS and strand for all DHSs (the 2 DHSs with the most reads removed); PCA of log(1 + count), centered and scaled; PC1 vs PC2 | `DHS_read_counts_PCA/` |

## Requirements

R: tidyverse, data.table, Rsamtools, reshape2, pbapply.

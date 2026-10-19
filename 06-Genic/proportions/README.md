# 06 — Genic: category proportions

Frequency of each 3′ end category per gene-associated DHS and genotype (`TLDRSeq_sup`), from
read counts scaled by the chrM scale factor of each replicate and pooled over the replicates of
each genotype. Input of the genic heatmap
([`../read-classes/03_read_categories_heatmap.r`](../read-classes)).

## Run

```bash
Rscript 02_build_freq_table.R
```

Run from this folder, after `../read-classes/label_reads_TLDRSeq_sup.R`. Paths are set in
[`../../config/paths.R`](../../config/paths.R).

## Inputs

| Input | From |
|---|---|
| reads per DHS, replicate and 3′ end label combination (`TTS_PROP_CATEGORIES_RDS`) | `../read-classes/label_reads_TLDRSeq_sup.R` |
| chrM read depth per replicate (`CHRM_DEPTH`) | 04-Normalization |

## Files

| File | Contents |
|---|---|
| `00_config.R` | input and output paths, the 21 replicates |
| `01_genic_proportion_curation.R` | `curate_proportion_categories()`: label combination → category + pA+ / pA− (pA+: poly(A) tail > 15 nt); `build_gene_count_long()`: counts per DHS, replicate and category, upstream reads and ZC3H4 removed; `calculate_pooled_counts_then_proportion()`: count / scale factor, summed over the replicates of each genotype, frequency |
| `02_build_freq_table.R` | scale factors (`define_scale_factor()`, [`04-Normalization/scale_factors.R`](../../04-Normalization)), then the curation and pooling of `01_genic_proportion_curation.R`; saves the table |

## Output

`GENIC_FREQ_RDS` (`$RESULTS_DIR/06-Genic/TLDRSeq_sup/proportions/`): one row per DHS × category × genotype
(WT, ARS2, RBM7, ZC3H18, ZCCHC8).

| Column | Content |
|---|---|
| `DHS` | DHS id |
| `feature_id` | category and pA status (`pAp` = pA+, `pAn` = pA−), or `others` |
| `geno` | genotype |
| `feature_count` | scaled read count of the category, summed over replicates |
| `geno_count` | scaled read count of all categories of the DHS in the genotype |
| `scaled_freq` | `feature_count / geno_count` |

## Requirements

R: tidyverse.

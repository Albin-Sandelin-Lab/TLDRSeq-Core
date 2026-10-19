# 04 — Normalization

TLDR-Seq read counts (`TLDRSeq_sup`) are normalized with one scale factor per replicate: the
replicate's chrM read depth divided by the mean chrM read depth over the 21 replicates.
Counts are divided by the scale factor.

## Run

```bash
Rscript 01_write_scale_factors.R
```

Paths are set in [`../config/paths.R`](../config/paths.R).

## Files

| File | Contents |
|---|---|
| `scale_factors.R` | `define_scale_factor()`, `validate_scale_factor_table()`, `extract_genotype_from_sample()`; the only definition of the scale factors, also used by `lib/deseq.R` (`run_deseq()`) and `06-Genic/proportions/` |
| `01_write_scale_factors.R` | writes the scale factor table |

## Input

chrM read depth per replicate (`CHRM_DEPTH`): one line per replicate, `<depth> <path/to/bed12_file>`.

## Output

`SCALE_FACTORS_CSV` (`$RESULTS_DIR/04-Normalization/scale_factor_per_rep.csv`): `sample`
(`scale_<replicate>`), `scale_factor`. Read by steps 05–07.

## Requirements

R: tidyverse.

#!/bin/bash
# PROMPT analyses (TLDRSeq_sup: WT, RBM7, ZCCHC8).
#
# Usage: ./run_all.sh
# Paths are set in ../config/paths.R

cd "$(dirname "${BASH_SOURCE[0]}")"

Rscript 01_prompt_reads_and_stats.r
Rscript 02_prompt_sensitivity.r
Rscript 03_prompt_vs_pmt_length.r
Rscript 04_prompt_length_categories.r
Rscript 05_prompt_pA_motifs.r
Rscript 06_prompt_windows.r
Rscript 07_prompt_AT_model.r
Rscript 08_prompt_gene_coverage.r

exit 0

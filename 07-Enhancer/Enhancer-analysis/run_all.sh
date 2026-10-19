#!/bin/bash
# Enhancer DHS analyses.
#
# Usage: ./run_all.sh
# Paths are set in ../../config/paths.R

cd "$(dirname "${BASH_SOURCE[0]}")"

Rscript 01_enhancer_coverage_and_chip.r
Rscript 02_enhancer_chip_logFC.r
for dataset in GSE142312 GSE115774; do
    Rscript 03_nascent_tss_edges.r $dataset
done
Rscript 04_enhancer_reads_sensitivity.r
Rscript 05_enhancer_tss_logo.r

exit 0

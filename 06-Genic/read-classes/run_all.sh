#!/bin/bash
# Genic read analyses.
#
# Usage: ./run_all.sh [n_cores]
# Paths are set in ../../config/paths.R

cd "$(dirname "${BASH_SOURCE[0]}")"

N_CORES=${1:-8}

## TLDRSeq_sup: 3' end labels, then category frequencies pooled per genotype
Rscript label_reads_TLDRSeq_sup.R $N_CORES
(cd ../proportions && Rscript 02_build_freq_table.R)

## nuclear / cytoplasmic and PAXT libraries
for dataset in TLDRSeq_nuclear_sup TLDRSeq_Z1_MTR4_sup; do
    Rscript 01_classify_reads_per_DHS.R $dataset $N_CORES
    Rscript 02_library_size_factors.r $dataset
done

Rscript 03_read_categories_heatmap.r
Rscript 04_cyto_nuc_ranks.r
Rscript 05_read_categories_sensitivity.r
Rscript 06_AT_model_polynomial.r

exit 0

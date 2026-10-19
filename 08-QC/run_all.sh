#!/bin/bash
# Quality control of the TLDRSeq_sup libraries (WT, RBM7, ZCCHC8).
#
# Usage: ./run_all.sh [n_cores]
# Paths are set in ../config/paths.R

cd "$(dirname "${BASH_SOURCE[0]}")"

N_CORES=${1:-15}

Rscript 01_library_size.r
Rscript 02_DHS_read_counts_PCA.r $N_CORES

exit 0

#!/bin/bash
# Define DHSs from the ENCODE mES DNase-seq data.
#
# Usage: ./run_all.sh ../config/DNase_ENCSR000CMW.sh

config=$(realpath "$1")
source "$config"
cd "$(dirname "${BASH_SOURCE[0]}")"

./01_download_dnase.sh "$config"
./02_map_dnase.sh "$config"
./03_call_peaks_dnase.sh "$config"
Rscript 04_decompose_DHS.r $DNASE_DIR $DHS_DIR

exit 0

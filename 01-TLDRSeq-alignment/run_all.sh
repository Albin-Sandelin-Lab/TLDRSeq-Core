#!/bin/bash
# Run the whole TLDR-Seq alignment step for one dataset.
#
# Usage: ./run_all.sh <dataset_config>
#   e.g. ./run_all.sh ../config/TLDRSeq_sup.sh
#        ./run_all.sh ../config/TLDRSeq_nuclear_sup.sh
#        ./run_all.sh ../config/TLDRSeq_Z1_MTR4_sup.sh

config=$(realpath "$1")
source "$config"
cd "$(dirname "${BASH_SOURCE[0]}")"

for library in "${LIBRARIES[@]}"
do
    ./01_trim_and_map.sh "$config" $library 12
    ./02_split_by_barcode.sh "$config" $library
done

./03_merge_per_sample.sh "$config"
./04_separate_pA_reads.sh "$config"

exit 0

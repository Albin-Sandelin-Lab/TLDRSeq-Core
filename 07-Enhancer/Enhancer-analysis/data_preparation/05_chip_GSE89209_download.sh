#!/bin/bash
# Download the GSE89209 ChIP-seq runs (mESC and AntNPC: H3K27me3, H3K27ac, p300, input) from SRA
# and merge the runs of each sample.
#
# Usage (from 07-Enhancer/Enhancer-analysis/data_preparation/): ./05_chip_GSE89209_download.sh
# Requires: sratoolkit (fasterq-dump), Entrez Direct (esearch, efetch)
# Output: $CHIP_GSE89209_FASTQ/<GSM>_<sample>[_1|_2].fastq

source ../../../config/paths.sh
module load sratoolkit

DATASETS=(
    "GSM2360930 WT_mESC_ChIP_H3K27me3"
    "GSM2416830 WT_AntNPC_ChIP_H3K27me3"
    "GSM2360929 WT_mESC_ChIP_H3K27ac"
    "GSM2360927 WT_AntNPC_ChIP_H3K27ac"
    "GSM2360934 WT_mESC_ChIP_p300"
    "GSM2433947 WT_AntNPC_ChIP_p300"
    "GSM2360933 WT_mESC_ChIP_input"
    "GSM2360928 WT_AntNPC_ChIP_input"
)

out_dir=$CHIP_GSE89209_FASTQ
mkdir -p "$out_dir"

for item in "${DATASETS[@]}"; do
    read -r GSM_ID NAME <<< "$item"
    echo "Processing $GSM_ID ($NAME)"

    # all SRR runs of this GSM
    SRR_LIST=$(esearch -db sra -query "${GSM_ID}" | efetch -format runinfo | cut -d ',' -f 1 | grep "^SRR")
    if [ -z "$SRR_LIST" ]; then
        echo "Error: No SRR files found for $GSM_ID!"
        continue
    fi
    echo "Found SRR runs for $GSM_ID: $SRR_LIST"

    r1_files=()
    r2_files=()
    single_files=()

    for SRR in $SRR_LIST; do
        echo "  Downloading $SRR..."
        fasterq-dump "$SRR" -e 10 --split-files -O "$out_dir"
        # paired-end or single-end
        if [[ -f "$out_dir/${SRR}_1.fastq" ]]; then
            r1_files+=("$out_dir/${SRR}_1.fastq")
            r2_files+=("$out_dir/${SRR}_2.fastq")
        elif [[ -f "$out_dir/${SRR}.fastq" ]]; then
            single_files+=("$out_dir/${SRR}.fastq")
        fi
    done

    # merge the runs into one FASTQ per sample
    if [ ${#r1_files[@]} -gt 0 ]; then
        cat "${r1_files[@]}" > "$out_dir/${GSM_ID}_${NAME}_1.fastq"
        cat "${r2_files[@]}" > "$out_dir/${GSM_ID}_${NAME}_2.fastq"
        rm "${r1_files[@]}" "${r2_files[@]}"
    elif [ ${#single_files[@]} -gt 0 ]; then
        cat "${single_files[@]}" > "$out_dir/${GSM_ID}_${NAME}.fastq"
        rm "${single_files[@]}"
    fi
done

exit 0

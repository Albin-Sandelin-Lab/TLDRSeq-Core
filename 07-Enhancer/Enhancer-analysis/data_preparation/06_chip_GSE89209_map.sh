#!/bin/bash
# Map the GSE89209 ChIP-seq samples (bowtie2, mm39), keep mapped reads and make bigWigs of the
# read coverage (deepTools bamCoverage, 10 bp bins, reads extended to fragments, no normalization:
# tracks are normalized by their area in 02_enhancer_chip_logFC.r).
#
# Usage (from 07-Enhancer/Enhancer-analysis/data_preparation/): ./06_chip_GSE89209_map.sh
# Requires: bowtie2, samtools, deeptools
# Input : $CHIP_GSE89209_FASTQ (05)
# Output: $CHIP_GSE89209_DIR/<sample>_filtered.bam(.bai), <sample>.bw (not for inputs)

source ../../../config/paths.sh
module load bowtie2
module load samtools
module load deeptools

data=$CHIP_GSE89209_FASTQ
results=$CHIP_GSE89209_DIR
mkdir -p "$results"

SAMPLES=(
    "GSM2360930_WT_mESC_ChIP_H3K27me3"
    "GSM2416830_WT_AntNPC_ChIP_H3K27me3"
    "GSM2360929_WT_mESC_ChIP_H3K27ac"
    "GSM2360927_WT_AntNPC_ChIP_H3K27ac"
    "GSM2360934_WT_mESC_ChIP_p300"
    "GSM2433947_WT_AntNPC_ChIP_p300"
    "GSM2360933_WT_mESC_ChIP_input"
    "GSM2360928_WT_AntNPC_ChIP_input"
)

# 1. Map and sort (paired-end or single-end)
for NAME in "${SAMPLES[@]}"; do
    echo "Mapping $NAME..."
    if [[ -f "$data/${NAME}_1.fastq" && -f "$data/${NAME}_2.fastq" ]]; then
        bowtie2 -x $BOWTIE2_INDEX -1 "$data/${NAME}_1.fastq" -2 "$data/${NAME}_2.fastq" -p 20 | \
        samtools sort -@ 20 -o "$results/${NAME}_sorted.bam"
    elif [[ -f "$data/${NAME}.fastq" ]]; then
        bowtie2 -x $BOWTIE2_INDEX -U "$data/${NAME}.fastq" -p 20 | \
        samtools sort -@ 20 -o "$results/${NAME}_sorted.bam"
    else
        echo "  Error: Fastq files for $NAME not found!"
    fi
done

# 2. Keep mapped reads, index
for NAME in "${SAMPLES[@]}"; do
    samtools view -h -F 4 -b "$results/${NAME}_sorted.bam" > "$results/${NAME}_filtered.bam"
    samtools index -M "$results/${NAME}_filtered.bam"
done

# 3. bigWig of the read coverage (not for the inputs)
for NAME in "${SAMPLES[@]}"; do
    if [[ "$NAME" != *"input"* ]]; then
        bamCoverage -b "$results/${NAME}_filtered.bam" \
                    -o "$results/${NAME}.bw" \
                    -p 20 \
                    --binSize 10 \
                    --extendReads
    fi
done

exit 0

#!/bin/bash
# Map the DNase-seq reads of each replicate with bwa mem.
#
# Usage: ./02_map_dnase.sh ../config/DNase_ENCSR000CMW.sh
#
# Output (in $DNASE_DIR): rep1.bam, rep2.bam (sorted, indexed)
#
# Requires: bwa, samtools

source "$1"

module load bwa
module load samtools

for rep in rep1 rep2
do
    bwa mem -t 20 $BWA_INDEX $DNASE_DIR/$rep.fq | samtools view -bS - | samtools sort -@20 | tee $DNASE_DIR/$rep.bam | samtools index - $DNASE_DIR/$rep.bam.bai
done

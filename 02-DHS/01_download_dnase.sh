#!/bin/bash
# Download the ENCODE DNase-seq reads (ENCSR000CMW, mES) and pool them per replicate.
#
# Usage: ./01_download_dnase.sh ../config/DNase_ENCSR000CMW.sh
#
# Output (in $DNASE_DIR): rep1.fq, rep2.fq

source "$1"

mkdir -p $DNASE_DIR/fastq
cd $DNASE_DIR/fastq

for f in "${DNASE_REP1[@]}" "${DNASE_REP2[@]}"
do
    wget https://www.encodeproject.org/files/$f/@@download/$f.fastq.gz &
done
wait

zcat $(printf "%s.fastq.gz " "${DNASE_REP1[@]}") > $DNASE_DIR/rep1.fq
zcat $(printf "%s.fastq.gz " "${DNASE_REP2[@]}") > $DNASE_DIR/rep2.fq

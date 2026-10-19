#!/bin/bash
# Call DNase-seq peaks with MACS3 (both replicates pooled) and make the signal bigWig.
#
# Usage: ./03_call_peaks_dnase.sh ../config/DNase_ENCSR000CMW.sh
#
# Output (in $DNASE_DIR):
#   test_peaks.narrowPeak   DNase peaks
#   test_treat_pileup.bw    DNase pileup signal
#
# Requires: MACS3 (3.0.0), bedGraphToBigWig

source "$1"

cd $DNASE_DIR

macs3 callpeak -t rep1.bam rep2.bam -g mm -n test -B -q 0.05

# Keep chr1-19, chrX and chrY only: the unplaced contigs are named differently in
# CHROM_SIZES (UCSC names) and would make bedGraphToBigWig fail; chrMT is removed too.
awk '$1 ~ /^chr([0-9]+|X|Y)$/' test_treat_pileup.bdg > test_treat_pileup_filtered.bdg

bedGraphToBigWig test_treat_pileup_filtered.bdg $CHROM_SIZES test_treat_pileup.bw

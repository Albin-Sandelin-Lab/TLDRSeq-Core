#!/bin/bash
# DNase-seq, mouse ES cells, ENCODE ENCSR000CMW (used to define DHSs)
source "$(dirname "${BASH_SOURCE[0]}")/paths.sh"

DNASE_DIR=$RESULTS_DIR/DNase_ENCSR000CMW

# ENCODE FASTQ files of each replicate
DNASE_REP1=(ENCFF001PIR ENCFF001PIS ENCFF001PIT ENCFF001PIU ENCFF001PIV)
DNASE_REP2=(ENCFF001PIW ENCFF001PIX)

# DHS definition (02-DHS/04_decompose_DHS.r). Folder name kept from the original analysis.
DHS_DIR=$RESULTS_DIR/dhs_charac/decompose_DHS

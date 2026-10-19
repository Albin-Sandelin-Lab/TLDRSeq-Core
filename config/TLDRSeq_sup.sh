#!/bin/bash
# TLDR-Seq, mES: WT, RBM7, ZCCHC8, ARS2, ZC3H4, ZC3H18 (see metadata.tsv)
source "$(dirname "${BASH_SOURCE[0]}")/paths.sh"

DATASET=TLDRSeq_sup
RAW_DIR=$DATA_DIR/TLDR-Seq_sup
METADATA=$RAW_DIR/metadata.tsv

# Number of barcodes used in the libraries (bc1..bcN)
N_BARCODES=21

# Basecalled (dorado, sup model) unaligned BAMs, one per sequencing library.
# lib5 is excluded (lib5_old.bam also exists); lib15_remade replaces lib15.
LIBRARIES=(
    basecall_lib1 basecall_lib2 basecall_lib3 basecall_lib4
    basecall_lib6 basecall_lib7 basecall_lib8 basecall_lib9 basecall_lib10
    basecall_lib11 basecall_lib12 basecall_lib13 basecall_lib14 basecall_lib15_remade
    basecall_lib16 basecall_lib17 basecall_lib18 basecall_lib19 basecall_lib20
    basecall_lib21
)

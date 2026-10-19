#!/bin/bash
# TLDR-Seq, mES: WT, MTR4, Z1-EV, Z1-FL, Z1-V712 (see metadata.tsv)
source "$(dirname "${BASH_SOURCE[0]}")/paths.sh"

DATASET=TLDRSeq_Z1_MTR4_sup
RAW_DIR=$DATA_DIR/TLDR-Seq_Z1_MTR4_sup
METADATA=$RAW_DIR/metadata.tsv

# Number of barcodes used in the library (bc1..bcN)
N_BARCODES=15

# Sample_Name in the metadata does not include the replicate (mES_WT = bc1, bc6, bc11):
# name the per-replicate BAMs <Sample_Name>_<Replicate> (e.g. mES_WT_rep1)
SAMPLE_NAME_ADD_REPLICATE=1

LIBRARIES=(basecall_lib1)

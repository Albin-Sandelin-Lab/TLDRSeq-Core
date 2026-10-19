#!/bin/bash
# TLDR-Seq, mES full / cytoplasmic / nuclear fractions
source "$(dirname "${BASH_SOURCE[0]}")/paths.sh"

DATASET=TLDRSeq_nuclear_sup
RAW_DIR=$DATA_DIR/TLDR-Seq_nuclear_sup
METADATA=$RAW_DIR/metadata.tsv

# Number of barcodes trimmed (bc1..bcN). The metadata only lists bc1-bc8;
# 9 is kept because it is what the published data were processed with.
N_BARCODES=9

LIBRARIES=(basecall_lib1)

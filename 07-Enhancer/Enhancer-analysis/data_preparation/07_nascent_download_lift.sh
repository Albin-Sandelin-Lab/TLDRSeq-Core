#!/bin/bash
# Download the strand-specific nascent RNA bigWigs (mm9) from GEO and lift them over to mm39
# with CrossMap.
#
# Usage (from 07-Enhancer/Enhancer-analysis/data_preparation/): ./07_nascent_download_lift.sh <GSE>
#   GSE142312 : GRO-seq, ESC and NPC (2 lines x 3 experiments)  -> $DATA_DIR/GSE142312_nascent/mm39/
#   GSE115774 : 4sU-seq, Day0 / 3 / 5 / 7 (2 replicates)        -> $DATA_DIR/GSE115774_nascent/bw39/
# Output files: <original name>.mm39.bw.bw (GSE142312), <original name>.bw (GSE115774)

source ../../../config/paths.sh
source ~/.bashrc
conda activate $CROSSMAP_ENV

GSE=$1
case $GSE in
    GSE142312) out_sub=mm39 ;;
    GSE115774) out_sub=bw39 ;;
    *) echo "GSE must be GSE142312 or GSE115774"; exit 1 ;;
esac

dir=$DATA_DIR/${GSE}_nascent
mkdir -p $dir
cd $dir

# GEO supplementary files: one bigWig per sample and strand
wget -O ${GSE}_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=${GSE}&format=file"
tar -xf ${GSE}_RAW.tar

# mm9 -> mm39
wget $CHAIN_MM9_MM39_URL
gzip -d mm9ToMm39.over.chain.gz

mkdir -p $out_sub
if [ "$GSE" == "GSE142312" ]; then
    ls | grep bw$ | xargs -I {} --max-proc=30 CrossMap bigwig mm9ToMm39.over.chain {} $out_sub/{}.mm39.bw
else
    ls | grep bigWig$ | xargs -I {} --max-proc=30 CrossMap bigwig mm9ToMm39.over.chain {} $out_sub/{}
fi

exit 0

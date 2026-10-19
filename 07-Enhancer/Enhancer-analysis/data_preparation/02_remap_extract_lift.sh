#!/bin/bash
# Extract the ReMap archive (bigwig_files/, mm10) and lift every bigWig over to mm39 with CrossMap.
#
# Usage (from 07-Enhancer/Enhancer-analysis/data_preparation/): ./02_remap_extract_lift.sh
# Input : $REMAP_DIR/ReMap_mmus_mm10_mESC_archive_part.tar.part.a{a..g} (01)
# Output: $REMAP_DIR/bigwig_files/<GSE>.<TF>.<biosample>[_<condition>].bigwig      (mm10)
#         $REMAP_DIR/bigwig_files_mm39/<same name>.bw                               (mm39)

source ../../../config/paths.sh
source ~/.bashrc
conda activate $CROSSMAP_ENV

cd $REMAP_DIR

# the 7 parts are consecutive pieces of one tar archive
cat ReMap_mmus_mm10_mESC_archive_part.tar.part.a? | tar -x

mkdir -p bigwig_files_mm39
ls bigwig_files | xargs -I {} --max-procs=20 CrossMap bigwig $CHAIN_MM10_MM39 bigwig_files/{} bigwig_files_mm39/{}

exit 0

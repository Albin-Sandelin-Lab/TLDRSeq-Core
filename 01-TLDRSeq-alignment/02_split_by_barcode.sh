#!/bin/bash
# Orient reads to the RNA strand and split one mapped library into one BAM per barcode.
#
# Usage: ./02_split_by_barcode.sh <dataset_config> <library>
#
# Input : $RESULTS_DIR/$DATASET/pipeline_map/<library>/mapping/trimmed_primary.bam
# Output: $RESULTS_DIR/$DATASET/get_primary_processed_bam/<library>/bcN.bam
#
# Read names carry the orientation and barcode added by cutadapt (s1_bcN_<read_id>,
# s2_bcN_<read_id>, see 01_trim_and_map.sh):
#   - only primary alignments are kept (flag 0 or 16)
#   - the strand of s2 (antisense) reads is flipped so that every read is on the RNA strand
#   - a read assigned to several barcodes is kept only once, in the barcode where its
#     MAPQ is the highest
#   - the s1/s2 and bcN prefixes are removed from the read name
#
# Requires: samtools

source "$1"
library=$2

module load samtools

data=$RESULTS_DIR/$DATASET/pipeline_map/$library/mapping
results=$RESULTS_DIR/$DATASET/get_primary_processed_bam/$library
mkdir -p $results

bam_genome=$data/trimmed_primary.bam
bam_header=$results/bam.header

samtools view -H $bam_genome > $bam_header

# After splitting the read name on "_": $1 = s1/s2, $2 = bcN, $3 = read id, $4 = flag, $7 = MAPQ
samtools view $bam_genome | awk -v FS="\t" -v OFS="\t" '{gsub("_","\t",$1)} ($2==0 || $2==16) {print $0}' |\
    awk -v OFS="\t" -v FS="\t" '{a=($1=="s2"?-1:1); b=($4==16?-1:1); $4=(a*b==1?0:16)} {print $0}' |\
    awk -v dir=$results -v OFS="\t" '(!($3 in max) || ($7 > max[$3])) {max[$3] = $7; fn[$3]=$2; $1=""; $2=""; line[$3] = substr($0,3)} END{for (i in line) {file=fn[i]".sam"; print line[i] > dir"/"file}}'

# One sorted and indexed BAM per barcode
find $results -name "*.sam" | xargs -I {} --max-procs=$N_BARCODES bash -c 'a=$(basename {} .sam); b=$(dirname {})/$a".bam"; cat "'$bam_header'" {} | samtools view -Sb | samtools sort | tee $b | samtools index - $b".bai"'

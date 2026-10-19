#!/bin/bash
# Merge the per-barcode BAMs of all libraries into one BAM per sample.
#
# Usage: ./03_merge_per_sample.sh <dataset_config>
#
# Input : $RESULTS_DIR/$DATASET/get_primary_processed_bam/<library>/bcN.bam
# Output: $RESULTS_DIR/$DATASET/concatenate_per_replicate/<Sample_Name>.bam
#
# Barcodes are matched to samples with the metadata (column 1 = barcode,
# column 7 = Sample_Name, column 8 = Replicate). For each barcode, the bcN.bam of all
# libraries are merged. Replicates are never merged together.
#
# The output is named <Sample_Name>, or <Sample_Name>_<Replicate> when the dataset
# config sets SAMPLE_NAME_ADD_REPLICATE=1 (datasets whose Sample_Name does not
# include the replicate, e.g. TLDRSeq_Z1_MTR4_sup).
#
# Requires: samtools

source "$1"

module load samtools

data=$RESULTS_DIR/$DATASET/get_primary_processed_bam
results=$RESULTS_DIR/$DATASET/concatenate_per_replicate
mkdir -p $results

# barcode <tab> output name
sample_list=$results/barcode_to_sample.tsv
awk -v OFS="\t" -v add_rep=${SAMPLE_NAME_ADD_REPLICATE:-0} 'NR>1 {print $1, (add_rep==1 ? $7"_"$8 : $7)}' $METADATA | sort -u > $sample_list

# Each output name must correspond to a single barcode, otherwise parallel jobs overwrite each other
if [ $(cut -f2 $sample_list | sort | uniq -d | wc -l) -ne 0 ]
then
    echo "Error: several barcodes have the same sample name in $METADATA:" >&2
    cut -f2 $sample_list | sort | uniq -d >&2
    exit 1
fi

cat $sample_list | xargs -I {} --max-procs=10 bash -c 'a=$(echo {} | awk "{print \$1}"); b=$(echo {} | awk "{print \$2}"); c=$(ls "'$data'"/*/$a.bam | tr "[\n]" "[ ]";); samtools merge -o /dev/stdout $c | samtools sort -@ 5 | tee "'$results'"/$b.bam | samtools index - "'$results'"/$b.bam.bai'

#!/bin/bash
# Split each sample BAM into reads with and without a poly(A) tail.
#
# Usage: ./04_separate_pA_reads.sh <dataset_config>
#
# Input : - basecalled BAMs $RAW_DIR/<library>.bam (poly(A) tail length estimated by
#           dorado, tag pt:i)
#         - $RESULTS_DIR/$DATASET/concatenate_per_replicate/<sample>.bam
# Output: $RESULTS_DIR/$DATASET/separate_pA_reads/<sample>_pA_plus.bam   (tail > 15 nt)
#         $RESULTS_DIR/$DATASET/separate_pA_reads/<sample>_pA_minus.bam  (other reads)
#
# Requires: samtools

source "$1"

module load samtools

data=$RESULTS_DIR/$DATASET/concatenate_per_replicate
results=$RESULTS_DIR/$DATASET/separate_pA_reads
mkdir -p $results

polyA_list=$results/reads_pA.txt
no_polyA_list=$results/no_polyA.list

# Poly(A) tail length of every read, from all libraries of the dataset
for library in "${LIBRARIES[@]}"
do
    samtools view $RAW_DIR/$library.bam | awk -v OFS="\t" '$1 !~ /^@/ {for (i=12; i<=NF; i++) if ($i ~ /^pt:i:/) print $1, $i}'
done > $results/reads_pA_length.tsv

# Reads with a tail > 15 nt, and all the other reads
sed 's/\tpt:i:/\t/' $results/reads_pA_length.tsv | awk '$2 > 15 {print $1}' > $polyA_list
for library in "${LIBRARIES[@]}"
do
    samtools view $RAW_DIR/$library.bam | awk '{print $1}'
done | cat $polyA_list - | awk '{count[$0]++} END {for (line in count) if (count[line] == 1) print line}' > $no_polyA_list

ls $data/*.bam | xargs -I {} --max-procs=10 bash -c 'b=$(basename {} .bam); c="'$results'"/${b}_pA_plus.bam; samtools view -h -b -N "'$polyA_list'" {} | samtools sort -@ 20 | tee $c | samtools index - ${c}.bai' &
ls $data/*.bam | xargs -I {} --max-procs=10 bash -c 'b=$(basename {} .bam); c="'$results'"/${b}_pA_minus.bam; samtools view -h -b -N "'$no_polyA_list'" {} | samtools sort -@ 20 | tee $c | samtools index - ${c}.bai'
wait

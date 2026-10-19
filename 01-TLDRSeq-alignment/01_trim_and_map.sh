#!/bin/bash
# Trim TLDR-Seq adapters/barcodes and map one sequencing library to the genome.
#
# Usage: ./01_trim_and_map.sh <dataset_config> <library> [steps]
#   dataset_config : e.g. ../config/TLDRSeq_sup.sh
#   library        : basename of the basecalled BAM in $RAW_DIR (e.g. basecall_lib1)
#   steps          : any combination of 1 (trim), 2 (map), 3 (one-sided trimming
#                    statistics, QC only). Default: 12
#
# Output (in $RESULTS_DIR/$DATASET/pipeline_map/<library>/):
#   all_pass.fastq.gz             reads extracted from the basecalled BAM
#   trimming/bcN_s{1,2}.fastq.gz  trimmed reads per barcode and orientation
#   trimmed.fastq.gz              all trimmed reads (read names prefixed s1_bcN_ / s2_bcN_)
#   mapping/trimmed_primary.bam   one alignment per read (highest AS)
#
# Requires: samtools, cutadapt, minimap2, GNU parallel, pigz

source "$1"
library=$2
steps=${3:-12}

module load samtools
module load cutadapt
module load minimap2

bam_in=$RAW_DIR/$library.bam
results=$RESULTS_DIR/$DATASET/pipeline_map/$library
results_trimmed=$results/trimming
results_trimmed_onesided=$results/trimming_onesided
results_mapped=$results/mapping
mkdir -p $results $results_trimmed $results_trimmed_onesided $results_mapped

fastq_in=$results/all_pass.fastq.gz

echo "Input: $bam_in"
echo "Results: $results"

##### Extract reads from the basecalled BAM
if [ ! -f $fastq_in ]
then
    samtools bam2fq -@ 40 $bam_in | pigz > $fastq_in
fi


##### Adapters / barcodes
# A2 adapter + barcode, one per barcode (bc1..bc21)
A2_bc_all=("TGGCTTGATCCCTCATCTTGTGAAGTTGTTTCGGTCGATTCCGTTTGTAGTCGTCTGT"
           "ACTGGTGCAGCTTTGAACATCTAGAGAGGGTACTATGTGCCTCAGCACCGTACAGCAA"
           "ATGGACTTTGGTAACTTCCTGCGTCACCCACACTTACTTCAGGACGTACCAGTAGAAG"
           "GTTGAATGAGCCTACTGGGTCCTCTTCTGAAGTTCCTGGGTCTTGAACCAGACTTGGT"
           "TGAGAGACAAGATTGTTCGTGGACGACAGACACCGTTCATCGACTTTCGGACGAAGAA"
           "AGATTCAGACCGTCTCATGCAAAGTTCTCAGTCTTCCTCCAGACAAGGCTTACGAAGC"
           "CAAGAGCTTTGACTAAGGAGCATGCCGATCCTTGTGGCTTCTAACTTCATGTCCCAGT"
           "TGGAAGATGAGACCCTGATCTACGGTTTGTCATACTCGTGTGCTCACCCTAGACACCT"
           "TCACTACTCAACAGGTGGCATGAAGAATCTAAGCAAACACGAAGGTGGTTGACAGACC"
           "GCTAGGTCAATCTCCTTCGGAAGTTACAGTCCGAGCCTCATGTGATCTACAACGTCAT"
           "CAGGTTACTCCTCCGTGAGTCTGAACCGAGATCCTACGAATGGAGTGTGTTGGGTAAC"
           "TCAATCAAGAAGGGAAAGCAAGGTCCTGGGAGCATCAGGTAGTAACAGGGAGGAAACA"
           "CATGTTCAACCAAGGCTTCTATGGTAGCTGACTGTCTTCCATACCGACAAGAAGGCAC"
           "TTGCGTCCTGTTACGAGAACTCATGCTTGCGATTGATGCTTAGTATCAGATCCAACAG"
           "GAGCCTCTCATTGTCCGTTCTCTAACCACAGGAGGACGATACAGAGAAGCTGTGTTCC"
           "ACCACTGCCATGTATCAAAGTACGCCACAGTGTCAACTAGAGCCTCTCGTGCAACTTT"
           "CTTACTACCCAGTGAACCTCCTCGTAGTTTGGATGACCAAGGATAGCCCATCTGGAAC"
           "GCATAGTTCTGCATGATGGGTTAGGGAGTTCGTCCAGAGAAGTACACGGCCTTCAGTG"
           "GTAAGTTGGGTATGCAACGCAATGCTACGTGTAAGGCATACCTGCCAGCATTCTCCTG"
           "CATACAGCGACTACGCATTCTCATCTTTCGTTGTTGACTCGACGGTAGCAGGTAGTTC"
           "CGACGGTTAGATTCACCTCTTACAAGTAGAAAGGGTTCCTTCCCACTCTACACCTGTA")

# Keep the barcodes used in this dataset
A2_bc=("${A2_bc_all[@]:0:$N_BARCODES}")
ix=($(seq 1 $N_BARCODES))

# A1 adapter (reverse complement)
A1rev=AAGAAAGTTGTCGGTGTCTTTGTGCGGATGAACATAGGATAGCGATTCACGGTATGTC

A2_bc_rev=($(printf "%s\n" "${A2_bc[@]}" | xargs -I {} bash -c "echo {} | tr '[ACGT]' '[TGCA]' | rev"))
A1=$(echo -n $A1rev | tr '[ACGT]' '[TGCA]' | rev)


##### Step 1: two-sided trimming, reads kept only if both adapters are found.
# s1 = read in sense orientation (A2+bc ... A1rev), s2 = antisense (A1 ... A2+bc rev)
if [[ $steps =~ "1" ]]
then
    echo "Trimming adapters (two-sided)"
    parallel -j $N_BARCODES --link "cutadapt -g {1}...{2} -o $results_trimmed/bc{3}_s1.fastq.gz -O 8 -m 1 --discard-untrimmed -x s1_bc{3}_ --buffer-size=16000000 $fastq_in > $results_trimmed/bc{3}_s1.log" ::: "${A2_bc[@]}" ::: "$A1rev" ::: "${ix[@]}" &
    parallel -j $N_BARCODES --link "cutadapt -g {1}...{2} -o $results_trimmed/bc{3}_s2.fastq.gz -O 8 -m 1 --discard-untrimmed -x s2_bc{3}_ --buffer-size=16000000 $fastq_in > $results_trimmed/bc{3}_s2.log" ::: "$A1" ::: "${A2_bc_rev[@]}" ::: "${ix[@]}"
    wait

    # Some cutadapt jobs occasionally fail (log with only 3 lines): rerun them until none fails
    nerr=$(wc -l $results_trimmed/*.log | awk '$1==3{print $2}' | wc -l)
    while [ ! $nerr -eq 0 ]
    do
        echo "Rerunning $nerr failed cutadapt jobs"
        cmd=$(wc -l $results_trimmed/*.log | awk '$1==3{print $2}' | parallel grep "parameter" {} | sed 's/Command line parameters:/cutadapt/' | awk '{print $0, "> "$5}' | sed 's/fastq.gz/log/3' | sed 's/-o/-j 50 -o/' | tr "[\n]" "[&]" | sed 's/&/ & /g')
        bash -c "$cmd wait"
        wait
        nerr=$(wc -l $results_trimmed/*.log | awk '$1==3{print $2}' | wc -l)
    done

    cat $results_trimmed/*.gz > $results/trimmed.fastq.gz
fi


##### Step 2: spliced mapping and selection of one alignment per read
if [[ $steps =~ "2" ]]
then
    echo "Mapping"
    bam=$results_mapped/trimmed.bam
    bam_H=$results_mapped/trimmed.headers
    bam_primary=$results_mapped/trimmed_primary.bam

    minimap2 -ax splice -un -t $THREADS $GENOME_FA $results/trimmed.fastq.gz | samtools view -S -b > $bam 2> $results_mapped/trimmed.log

    # Remove secondary (256) and unmapped (4) alignments, then keep for each read the
    # alignment with the highest alignment score (AS, column 14)
    samtools view -H $bam > $bam_H
    samtools view $bam -F 256 -F 4 | awk -v OFS="\t" '{gsub("AS:i:","",$14)}{print $0}' |\
        awk -v OFS="\t" '(!($1 in max) || ($14 > max[$1])) {max[$1] = $14; line[$1] = $0} END{for(i in line) print line[i]}' |\
        awk -v OFS="\t" '{gsub("^","AS:i:",$14)} {print $0}' |\
        cat $bam_H - | samtools view -S -b | samtools sort -@ 45 | tee $bam_primary | samtools index - $bam_primary.bai
fi


##### Step 3 (QC only): one-sided trimming, used only for the cutadapt statistics in the logs
if [[ $steps =~ "3" ]]
then
    echo "Trimming adapters (one-sided, statistics only)"
    parallel -j $N_BARCODES --link "cutadapt -a {1}...{2} -O 8 -m 1 --discard-untrimmed -x s1_bc{3}_ --buffer-size=16000000 $fastq_in 2> $results_trimmed_onesided/bc{3}_s1.log | cat > /dev/null" ::: "${A2_bc[@]}" ::: "$A1rev" ::: "${ix[@]}" &
    parallel -j $N_BARCODES --link "cutadapt -a {1}...{2} -O 8 -m 1 --discard-untrimmed -x s2_bc{3}_ --buffer-size=16000000 $fastq_in 2> $results_trimmed_onesided/bc{3}_s2.log | cat > /dev/null" ::: "$A1" ::: "${A2_bc_rev[@]}" ::: "${ix[@]}"
    wait
fi

exit 0

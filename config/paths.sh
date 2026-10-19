#!/bin/bash
# Shared paths for the whole pipeline. Edit this file to run on another system.
# Sourced by every script (via a dataset config in this folder).

# Root folder of the raw data (basecalled unaligned BAMs, metadata, external data)
DATA_DIR=/maps/projects/sandelin_main/data/degradation_project

# Root folder where every step writes its results (<RESULTS_DIR>/<DATASET>/<step>/)
RESULTS_DIR=/maps/projects/sandelin_main/people/fzh976/degradation_project/results

# Mouse genome (GRCm39 / mm39 primary assembly, UCSC-style chromosome names: chr1, chr2, ...)
GENOME_FA=/maps/projects/sandelin_main/data/shared_data/genomes/mm39/Mus_musculus.GRCm39.dna.primary_assembly.fa

# Temporary files (samtools sort, etc.)
export TMPDIR=/maps/projects/sandelin_main/scratch

# Threads
THREADS=50

# bwa index of the same genome (DNase-seq mapping)
BWA_INDEX=/maps/projects/sandelin_main/data/shared_data/bwa_indexes/mm39/Mus_musculus.GRCm39.dna.primary_assembly.fa

# Chromosome sizes (bedGraph -> bigWig)
CHROM_SIZES=/maps/projects/sandelin_main/data/shared_data/genomes/mm39/mm39.chrom.sizes

#### External data (07-Enhancer/Enhancer-analysis/data_preparation)
# bowtie2 index of the genome (ChIP-seq mapping)
BOWTIE2_INDEX=/maps/projects/sandelin_main/data/shared_data/bowtie2_indexes/mm39/mm39
# liftOver chain files (CrossMap)
CHAIN_MM10_MM39=/maps/projects/sandelin_main/data/shared_data/chain_files/mm10ToMm39.over.chain.gz
CHAIN_MM9_MM39_URL=https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm39.over.chain.gz
# conda environment with CrossMap
CROSSMAP_ENV=crossmap
# ReMap 2022 mESC bigWigs (mm10), Recherche Data Gouv doi:10.57745/GEIHVN
REMAP_DIR=$DATA_DIR/bw_remap
# GSE89209 ChIP-seq FASTQs (mESC / AntNPC) and processed files
CHIP_GSE89209_FASTQ=$DATA_DIR/GSE89209_ChIP/fastq_files
CHIP_GSE89209_DIR=$RESULTS_DIR/diverse/ChIP_GSE89209

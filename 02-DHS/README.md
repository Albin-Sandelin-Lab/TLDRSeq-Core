# 02 — DHS definition

DNase I hypersensitive sites (DHSs) are defined from ENCODE DNase-seq of mouse ES cells
([ENCSR000CMW](https://www.encodeproject.org/experiments/ENCSR000CMW/)). MACS3 peaks
are merged, then split into individual DHSs at the valleys of the smoothed DNase signal.

## Run

```bash
./run_all.sh ../config/DNase_ENCSR000CMW.sh
```

## Steps

| Script | Does | Output |
|---|---|---|
| `01_download_dnase.sh` | download ENCODE FASTQs; pool rep1 = ENCFF001PIR, PIS, PIT, PIU, PIV and rep2 = ENCFF001PIW, PIX | `$DNASE_DIR/rep{1,2}.fq` |
| `02_map_dnase.sh` | `bwa mem` on GRCm39 | `$DNASE_DIR/rep{1,2}.bam` |
| `03_call_peaks_dnase.sh` | `macs3 callpeak -t rep1.bam rep2.bam -g mm -B -q 0.05`; pileup restricted to chr1–19, X, Y → bigWig | `$DNASE_DIR/test_peaks.narrowPeak`, `test_treat_pileup.bw` |
| `04_decompose_DHS.r` | smooth DNase signal (running mean 100 bp, 50 bp bins); extend peaks ±600 bp and merge; split merged regions at signal valleys (min distance 250 bp, summit > 0.3 × max, summit > 1.6 × (valley + 1)); score = DNase coverage | `$DHS_DIR/decomposed_peaks.bed`, `decomposed_peaks_score.bed`, `processed_peaks.bed`, `smoothed_dnase1.bw` |

`decomposed_peaks_score.bed` is the DHS list used in the following steps.

## Requirements

wget, bwa, samtools, MACS3 3.0.0, bedGraphToBigWig (UCSC), R with tidyverse,
GenomicRanges, rtracklayer, plyranges, BiocParallel, data.table.

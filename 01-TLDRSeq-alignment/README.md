# 01 — TLDR-Seq alignment

From basecalled nanopore reads to one BAM per sample, with reads oriented to the RNA
strand, and the same BAMs split into reads with and without a poly(A) tail.

The same scripts are used for the three TLDR-Seq datasets; the dataset is selected
by a config file in [`../config/`](../config):

| Config | Libraries | Barcodes |
|---|---|---|
| `TLDRSeq_sup.sh` | 20 | 21 |
| `TLDRSeq_nuclear_sup.sh` | 1 | 9 |
| `TLDRSeq_Z1_MTR4_sup.sh` | 1 | 15 |

Shared paths (data, results, genome, temp folder) are set in `../config/paths.sh`.

## Run

```bash
./run_all.sh ../config/TLDRSeq_sup.sh
```

## Steps

| Script | Runs on | Does | Output (in `$RESULTS_DIR/<dataset>/`) |
|---|---|---|---|
| `01_trim_and_map.sh` | each library | BAM → FASTQ; two-sided cutadapt trimming per barcode and orientation (s1 = sense, s2 = antisense); `minimap2 -ax splice -un` on GRCm39; one alignment per read (highest `AS`) | `pipeline_map/<lib>/mapping/trimmed_primary.bam` |
| `02_split_by_barcode.sh` | each library | keep flag 0/16, flip strand of s2 reads, keep each read in one barcode only (highest MAPQ), one BAM per barcode | `get_primary_processed_bam/<lib>/bcN.bam` |
| `03_merge_per_sample.sh` | dataset | merge each barcode across libraries (replicates are kept separate), named by `Sample_Name` (+ `_<Replicate>` if `SAMPLE_NAME_ADD_REPLICATE=1`) | `concatenate_per_replicate/<sample>.bam` |
| `04_separate_pA_reads.sh` | dataset | split samples on the dorado poly(A) length (`pt:i` tag, > 15 nt) | `separate_pA_reads/<sample>_pA_{plus,minus}.bam` |

`01_trim_and_map.sh` also has an optional step `3` (one-sided trimming), used only to
produce cutadapt statistics for QC.

## Requirements

samtools, cutadapt, minimap2, GNU parallel, pigz (loaded with `module load` on our cluster).

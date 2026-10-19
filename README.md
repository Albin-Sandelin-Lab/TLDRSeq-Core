# TLDRSeq-Core

Source code for the analysis in **Full-length transcript architecture reveals RNA fate and regulatory state**.

TLDR-seq (**T**rue **l**ength of **d**iverse capped **R**NA sequencing) is a full-length RNA sequencing method based on nanopore sequencing with cap trapping. It sequences polyadenylated and non-polyadenylated capped RNAs end-to-end, so each read can be classified by where its 3' end falls relative to gene structure.

The transcripts are studied at DNase I hypersensitive sites (DHSs) in mouse embryonic stem cells, together with their sensitivity to depletion of exosome cofactors (NEXT: RBM7, ZCCHC8; PAXT: MTR4, Z1). Genome: mouse GRCm39 (mm39).

## Repository

```
01-TLDRSeq-alignment/       nanopore reads -> one strand-oriented BAM per sample; pA+ / pA- reads
02-DHS/                     DHS calling, sample universe, BED12 reader
03-Annotation/              DHS-to-gene maps, RefSeq Select annotation
04-Normalization/           scale factors derived from mitochondrial read depth
05-PROMPT/                  promoter upstream transcripts: sensitivity, length, pA status, 3' end motifs, AT content
06-Genic/
  peak-calling/             the representative TTS of each gene
  read-classes/             how each read's 3' end is classified;
                            nuclear / cytoplasmic localization;
                            sensitivity per class; AT content model
  proportions/              class curation, scale-factor division, pooling
  clustering/               genes clustered on their class proportions
  covariates/               the sequence, conservation and codon covariates
  regression/               the CLR linear model
  gsea/                     GO-term GSEA over class proportions
07-Enhancer/
  class-definition/         enhancer DHS classification, enhancer read summaries
  erna-characterization/    eRNA length, TSS spacing, STARR-seq and TF density
  Enhancer-analysis/        enhancer coverage, chromatin and TF signal, nascent RNA, sensitivity,
                            TSS sequence; data_preparation/ for the public datasets
08-QC/                      library size, PCA of the read counts per DHS
config/                     paths (paths.sh, paths.R) and one config per sequencing dataset
lib/                        shared R helpers (DESeq2 with chrM size factors, reads per DHS,
                            peak calling, PWM scanning) and MEME -> PFM conversion
```

Folders with a `README.md` and a `run_all.sh` (`01-TLDRSeq-alignment`, `02-DHS`, `05-PROMPT`,
`06-Genic/read-classes`, `07-Enhancer/Enhancer-analysis`, `08-QC`) list their inputs, scripts, outputs and
requirements, and run their scripts in order.

## Datasets

| Config | Data |
|---|---|
| `config/TLDRSeq_sup.sh` | TLDR-Seq, WT and RBM7, ZCCHC8, ARS2, ZC3H18, ZC3H4 depleted mESCs |
| `config/TLDRSeq_nuclear_sup.sh` | TLDR-Seq, total, nuclear and cytoplasmic RNA |
| `config/TLDRSeq_Z1_MTR4_sup.sh` | TLDR-Seq, WT, MTR4 and Z1 lines |
| `config/DNase_ENCSR000CMW.sh` | ENCODE DNase-seq of mESCs |

## Running

1. Edit [`config/paths.sh`](config/paths.sh) and [`config/paths.R`](config/paths.R): data, results,
   genome and annotation paths.
2. Run the steps in order. Shell steps take a dataset config as first argument; R scripts are run
   from their own folder:

```bash
(cd 01-TLDRSeq-alignment && ./run_all.sh ../config/TLDRSeq_sup.sh)   # and the other TLDR-Seq configs
(cd 02-DHS && ./run_all.sh ../config/DNase_ENCSR000CMW.sh)
(cd 04-Normalization && Rscript 01_write_scale_factors.R)
(cd 05-PROMPT && ./run_all.sh)
(cd 06-Genic/read-classes && ./run_all.sh)
(cd 07-Enhancer/Enhancer-analysis && ./run_all.sh)                    # after data_preparation/
(cd 08-QC && ./run_all.sh)
```

`06-Genic/read-classes/06_AT_model_polynomial.r` uses outputs of `05-PROMPT`. Results are written
to `$RESULTS_DIR/<dataset>/…` (01, 02) and `$RESULTS_DIR/<step>/…` (05–08).

Conventions: replicates are kept separate (no pooling of BAM files); pA+ reads have a poly(A)
tail > 15 nt (dorado estimate); TLDR-Seq read counts are normalized with chrM-based scale factors.

## Requirements

Command line: dorado-basecalled BAMs as input, samtools, cutadapt, minimap2, GNU parallel, pigz,
bwa, MACS3, bedGraphToBigWig, MEME suite (XSTREME), bowtie2, deepTools, CrossMap, sratoolkit.

R (≥ 4.2): tidyverse, data.table, GenomicRanges, GenomicAlignments, GenomicFeatures, rtracklayer,
plyranges, Biostrings, BSgenome, DESeq2, apeglm, bambu, classInt (with e1071), ggseqlogo,
pheatmap, patchwork and the packages listed in each folder's README.

## Reference

Auxillos J, Stigliani A, Vaagensø CS, Garland W, Niazi AM, Valen E, Jensen TH, Sandelin A. True length of diverse capped RNA sequencing (TLDR-seq): 5'-3'-end sequencing of capped RNAs regardless of 3'-end status. *Nucleic Acids Res.* 2025; 53(6):gkaf240. doi:10.1093/nar/gkaf240

## Code availability

This repository contains the analysis source code accompanying the manuscript "Full-length transcript architecture reveals RNA fate and regulatory state" (Dai Y, Stigliani A, Auxillos JY, Garland W, Williams K, Souville L, Ballester B, Amundarain A, Parviz M, George A, Vaagensø C, Andersson R, Jensen TH, Sandelin A; submitted). 

## License

MIT. See `LICENSE`.

## Define DHSs from the DNase-seq peaks and signal.
##
## Usage: Rscript 04_decompose_DHS.r <dnase_dir> <out_dir> [n_cores]
##   dnase_dir : folder with test_peaks.narrowPeak and test_treat_pileup.bw (03_call_peaks_dnase.sh)
##   out_dir   : output folder
##
## 1. The DNase signal is smoothed (running mean over 100 bp, averaged in 50 bp bins).
## 2. MACS3 peaks are extended by 600 bp on each side and overlapping peaks are merged.
## 3. Each merged region is split into individual DHSs at the valleys of the smoothed
##    signal (find_local_maxima.3).
## 4. Each DHS is scored by the total smoothed DNase signal it covers.
##
## Output (in out_dir):
##   smoothed_dnase1.bw          smoothed DNase signal
##   processed_peaks.bed         merged extended peaks
##   decomposed_peaks.bed        DHSs
##   decomposed_peaks_score.bed  DHSs with a DNase signal, score = DNase coverage

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(plyranges)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
dnase_dir <- args[1]
out_dir <- args[2]
n_cores <- ifelse(length(args) >= 3, as.integer(args[3]), 75)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dhs.bw.name <- file.path(dnase_dir, "test_treat_pileup.bw")
dhs.bed.name <- file.path(dnase_dir, "test_peaks.narrowPeak")

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

dhs.bed <- import(dhs.bed.name, format = "BED", extraCols = extraCols_narrowPeak)
dhs.bw <- import.bw(dhs.bw.name)
seqlevels(dhs.bed, pruning.mode="coarse") <- seqlevels(dhs.bw)


#### 1. Smooth the DNase signal

bins <- tileGenome(seqinfo(dhs.bw), tilewidth = 50, cut.last.tile.in.chrom = TRUE)

## Running mean over k bp, then average score per bin
k <- 100
rle <- coverage(dhs.bw, weight = "score")
rle_mean <- BiocParallel::bplapply(rle, runmean, k = k,
                                   BPPARAM = BiocParallel::MulticoreParam(workers = 21))

gr_binned <- binnedAverage(bins, as(rle_mean, "RleList"), "score")

## runmean drops the first k-1 positions: shift by k/2 to re-center the signal
gr_binned.filtered <- gr_binned[gr_binned$score != 0] %>%
    mutate(end = end + floor(k/2), start = start + floor(k/2))

export.bw(gr_binned.filtered, file.path(out_dir, "smoothed_dnase1.bw"))


#### 2. Extend and merge the peaks

processed.peaks <- dhs.bed %>%
    mutate(start = start - 600, end = end + 600) %>%
    reduce()

export(processed.peaks, file.path(out_dir, "processed_peaks.bed"))


#### 3. Split the merged regions into DHSs

## vec          : smoothed signal at each position of the region
## min_distance : minimum distance between two DHS summits
## fraction     : summits must be above fraction * max(vec)
## coeff        : a valley separates two summits if summit > coeff * (valley + 1)
## start, chrom : genomic coordinates of the region
find_local_maxima.3 <- function(vec, min_distance, fraction, coeff, start = 0, peak_id, chrom) {
    ## Find global maximum and threshold
    max.vec <- max(vec)
    th <- fraction * max.vec
    ## Remove duplicates and calculate midpoints
    a.score <- vec[!duplicated(data.table::rleid(vec))]
    a.rep <- diff(which(!duplicated(data.table::rleid(c(vec, Inf)))))
    a.rep.cumul <- cumsum(a.rep)
    a.mid <- a.rep.cumul - floor(a.rep/2)
    ## Identify all local maxima above fraction * max.vec
    maxima_idx <- which((a.score > th) & diff(sign(diff(c(-Inf, a.score, -Inf)))) == -2)
    maxima_idx <- maxima_idx[!maxima_idx %in% c(1, a.score %>% length)] ## remove if maxima are at the edges of the peak
    summits <- a.mid[maxima_idx]
    summit_values <- vec[summits]
    ## Identify all local minima
    minima_idx <- which(diff(sign(diff(c(Inf, a.score, Inf)))) == 2)
    if (length(minima_idx) == 0) {
        minima <- c(1, length(vec))  ## Fallback to endpoints
    } else {
        minima <- a.mid[minima_idx]
    }
    minima_values <- vec[minima]
    ## Find position of global maximum
    global_max_idx <- summits[which.max(summit_values)]
    ## Sort minima by distance to global maximum
    minima_order <- order(abs(minima - global_max_idx))
    sorted_minima <- minima[minima_order]
    sorted_minima_values <- minima_values[minima_order]
    ## Filter summits based on minima
    final_summits <- integer(0)
    candidates.all <- integer(0)
    flagged.all <- integer(0)
    for (i in seq_along(sorted_minima)) {
        summits <- summits[!summits %in% flagged.all]
        m <- sorted_minima[i]
        m_val <- sorted_minima_values[i]
        min_threshold <- coeff * (m_val + 1)
        ## Find closest left and right summits
        left_summits <- summits[rev(which(summits < m))[1]]
        right_summits <- summits[which(summits > m)[1]]
        if(is.na(left_summits))
            left_summits <- NULL
        if(is.na(right_summits))
            right_summits <- NULL
        left_summit <- if (length(left_summits) > 0) left_summits[which.max(left_summits)] else NULL
        right_summit <- if (length(right_summits) > 0) right_summits[which.min(right_summits)] else NULL
        left_val <- if (!is.null(left_summit)) vec[left_summit] else -Inf
        right_val <- if (!is.null(right_summit)) vec[right_summit] else -Inf
        ## Evaluate summits against threshold
        left_ok <- !is.null(left_summit) && left_val > min_threshold
        right_ok <- !is.null(right_summit) && right_val > min_threshold
        candidates <- integer(0)
        if (left_ok && right_ok) {
            ## Keep both if above threshold
            candidates <- c(left_summit, right_summit)
        } else if (left_ok) {
            ## Keep left if only it's above
            candidates <- left_summit
            flagged.all <- c(flagged.all, right_summit)
        } else if (right_ok) {
            ## Keep right if only it's above
            candidates <- right_summit
            flagged.all <- c(flagged.all, left_summit)
        } else if ( !is.null(left_summit) || !is.null(right_summit)){
            if (left_val > right_val) {
                candidates <- left_summit
                flagged.all <- c(flagged.all, right_summit)
            } else {
                candidates <- right_summit
                flagged.all <- c(flagged.all, left_summit)
            }
        }
        candidates.all <- c(candidates.all, candidates)
    }
    candidates.all <- unique(candidates.all)
    flagged.all <- unique(flagged.all)
    candidates.all <- candidates.all[!candidates.all %in% flagged.all]
    ## Add candidates if they respect min_distance
    sorted_summits <- summits[order(-summit_values)]
    sorted_summits_filtered <- sorted_summits[sorted_summits %in% candidates.all]
    for (p in sorted_summits_filtered) {
        if (all(abs(final_summits - p) >= min_distance)) {
            final_summits <- c(final_summits, p)
        }
    }
    ## DHS boundaries = deepest minimum between two consecutive summits
    m1 <- c(1, minima, length(vec)) %>% unique
    m1.val <- vec[m1]
    s <- c(0, sort(final_summits), length(vec) +1) %>% unique
    start.peaks <- numeric()
    end.peaks <- numeric()
    for (i in 2 : (length(s) -1 ) ){
        m1.left <- m1[m1 > s[i - 1] & m1 < s[i]][which.min(m1.val[m1 > s[i - 1] & m1 < s[i]])]
        m1.right <- m1[m1 > s[i] & m1 < s[i + 1]][which.min(m1.val[m1 > s[i] & m1 < s[i + 1]])]
        if(i == 2){
            m1.left <- minima[1]
        }
        if (i == (length(s) -1)){
            m1.right <- minima[length(minima)]
        }
        start.peaks <- c(start.peaks, m1.left + start)
        end.peaks <- c(end.peaks, m1.right + start)
    }
    if(length(start.peaks) != length(end.peaks)){
        print(peak_id)
        stop("Error")
    }
    return(data.frame(seqnames = rep(chrom, length(end.peaks)), start = start.peaks, end = end.peaks))
}

dhs.bw <- gr_binned.filtered

## Smoothed signal at each position of each merged region
overlaps <- findOverlaps(processed.peaks, dhs.bw)
bw_subset <- dhs.bw[subjectHits(overlaps)]
regions_subset <- processed.peaks[queryHits(overlaps)]

bw_clipped <- pintersect(bw_subset, regions_subset, ignore.strand = TRUE)

all_positions <- unlist(tile(processed.peaks, width = 1))
all_positions$peak_id <- rep(1:length(processed.peaks), width(processed.peaks))

bw_clipped$peak_id <- queryHits(overlaps)
signal_hits <- findOverlaps(all_positions, bw_clipped, type = "within")
all_positions$score <- 0
all_positions$score[queryHits(signal_hits)] <- bw_clipped$score[subjectHits(signal_hits)]

all_position.score.list <- split(all_positions, all_positions$peak_id) %>%
    mclapply(function(x) mcols(x)$score, mc.cores = n_cores)
peak.starts <- processed.peaks %>% start()

decomposed_peaks <- mcmapply(find_local_maxima.3, all_position.score.list,
                             min_distance = 250, fraction = 0.3,
                             coeff = 1.6, start = peak.starts, peak_id = 1:length(peak.starts),
                             chrom = seqnames(processed.peaks) %>% as.character(),
                             mc.cores = n_cores,
                             SIMPLIFY = FALSE)

bed.processed.peaks <- decomposed_peaks %>%
    do.call(rbind.data.frame , .) %>%
    as_granges() %>%
    mutate(name = 1:length(.))
export(bed.processed.peaks, file.path(out_dir, "decomposed_peaks.bed"))


#### 4. DNase coverage over each DHS

overlaps2 <- findOverlaps(bed.processed.peaks, dhs.bw)
bw_subset.2 <- dhs.bw[subjectHits(overlaps2)]
regions_subset.2 <- bed.processed.peaks[queryHits(overlaps2)]

bw_clipped.2 <- pintersect(bw_subset.2, regions_subset.2, ignore.strand = TRUE)
bw_clipped.2$peak_id <- queryHits(overlaps2)

## DHSs without any signal are dropped
bed.processed.peaks_filtered <- bed.processed.peaks[unique(queryHits(overlaps2))]
bed.processed.peaks_filtered$peak_id <- unique(queryHits(overlaps2))

bw_clipped.2$coverage <- bw_clipped.2$score * width(bw_clipped.2)

bed.processed.peaks.cov <- aggregate(bw_clipped.2$coverage,
                                     by = list(bw_clipped.2$peak_id),
                                     FUN = sum) %>%
    dplyr::rename(peak_id = Group.1, score = x)

bed.processed.peaks_filtered$score <- bed.processed.peaks.cov$score

export(bed.processed.peaks_filtered, file.path(out_dir, "decomposed_peaks_score.bed"))

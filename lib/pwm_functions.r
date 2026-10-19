make.PWM <- function(x)
{
    nRegion <- rowSums(x[1,])
    x.t <- t(as.matrix(x))
    max_freq <- apply(x.t, FUN=max, 2)
    maxi_freq <- matrix(nrow=4, rep(max_freq, 4), byrow=TRUE)
    pwm <- log(x.t/maxi_freq)
    return(pwm)
}


get_pwm_score <- function(pwm, seq_chars, start, len) {
  nts <- seq_chars[start:(start + len - 1)]
  # Lookup table for nucleotide-to-index mapping
  nt_to_idx <- c(A = 1, C = 2, G = 3, T = 4)
  # Map nucleotides to indices
  indices <- nt_to_idx[nts]
  score <- sum(pwm[cbind(indices, seq_len(len))], na.rm = TRUE)  
  return(score)
}

get_max_pwm_score <- function(sequence, pwm, position.mode = FALSE, threshold = 0) {
  pwm_len <- ncol(pwm)
  seq_len <- nchar(sequence)
  stopifnot(seq_len >= pwm_len)  
  # Pre-split sequence into character vector once
  seq_chars <- strsplit(sequence, "")[[1]]
  # Slide window and compute scores
  n <- seq_len - pwm_len + 1
  max_score <- -Inf
  if(position.mode){
      score <- numeric(n)
      for (i in seq_len(n)) {
          score[i] <- get_pwm_score(pwm, seq_chars, i, pwm_len)
      }
      if(threshold == 0){
          return(which( score == max(score)))
      }else{
          min.pwm.score = apply(pwm, 2, min)
          threshold.value = min.pwm.score/(1/threshold)
          return(which(score >= threshold.value))
      }
  }else{
      for (i in seq_len(n)) {
          score <- get_pwm_score(pwm, seq_chars, i, pwm_len)
          if (score > max_score) max_score <- score
          if (max_score == 0) {
              break  # Early exit if score is zero
          }   
      }
      return(max_score)
  }
}


run_1_pwm_on_seq <- function(pwm, seqs, name.pwm, position.mode, threshold) {
    scores <- mclapply(seqs, get_max_pwm_score, pwm, position.mode, threshold, mc.cores = 15) %>%
        unlist()
    if(position.mode){
        return(scores)
    }
    scores_df <- data.frame(seqname = names(seqs), score = scores, ID = rep(name.pwm, length(seqs)))    
    return(scores_df)
}


run_all_pwm_on_seq <- function(pwm_list, seqs, position.mode, threshold) {
    name.pwm <- names(pwm_list)
    scores_list <- mcmapply( function(x, y, z, p, t) run_1_pwm_on_seq(x, y, z, p, t),
                            x = pwm_list,
                            z = name.pwm,
                            MoreArgs = list(y = seqs, p = position.mode, t = threshold),
                            SIMPLIFY=FALSE,
                            mc.cores = 5)
    return(scores_list)
}


compute_AUC <- function(df.pos, df.neg, name.elt) {
    scores.pos <- df.pos$absScore
    scores.neg <- df.neg$absScore
    w <- wilcox.test(scores.pos, scores.neg, alternative = "greater")
    pval <- w$p.value
    D <- w$statistic
    AUC <- D/ (length(scores.pos) * length(scores.neg))
    return(data.frame(AUC=AUC, pval = pval, name.elt =name.elt))
}

make.pos.neg <- function(list.df, name.pwm)
{
    names.clusters <- names(list.df)
    a <- lapply(names.clusters, function(x, n, l) {
        pos <- l[[x]]$score %>% 
            data.frame(absScore = .)
        neg <- l[names(l) != x] %>%
            lapply(getElement, "score") %>%
            unlist %>%
            data.frame(absScore = .)
        return(list(pos, neg, n))
    },
    name.pwm, list.df)
    names(a) <- names.clusters
    return(a)
}


make.pos.neg.1.1 <- function(list.df, name.pwm) {
    n <- names(list.df)
    l <- unlist(lapply(n, function(nx) {
        m <- n[n != nx]
        lapply(m, function(mx) {
            pos <- data.frame(absScore = list.df[[nx]]$score)
            neg <- data.frame(absScore = list.df[[mx]]$score)
            list(pos, neg, nx, mx, name.pwm)
        })
    }), recursive = FALSE)
    return(l)
}

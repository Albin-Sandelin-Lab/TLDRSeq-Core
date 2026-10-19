# exon_structure.R -- parse exon coordinates and locate positions in a transcript

# ---- parse_exon_structure --------------------------------------------------
# Parse comma-separated exon coordinates into a validated, structured list.
#
# Parameters:
#   exon_starts — character (comma-separated) or numeric vector of exon starts
#   exon_ends   — character (comma-separated) or numeric vector of exon ends
#   strand      — "+" or "-"
#   validate    — logical, perform validation checks (default TRUE)
#
# Returns:
#   list with:
#     exons             — data.frame (exon_num, genomic_idx, start, end, length, is_first, is_last)
#                         exon_num follows TRANSCRIPT order (1 = first in 5'->3')
#     strand            — character
#     gene_start        — leftmost genomic coordinate
#     gene_end          — rightmost genomic coordinate
#     total_exonic_length — sum of all exon lengths
#     n_exons           — integer
#     is_single_exon    — logical
parse_exon_structure <- function(
  exon_starts,
  exon_ends,
  strand,
  validate = TRUE
) {
  # Handle NA/NULL
  if (is.null(exon_starts) || is.null(exon_ends) || is.null(strand)) {
    return(NULL)
  }
  if (is.na(exon_starts) || is.na(exon_ends) || is.na(strand)) {
    return(NULL)
  }

  # Parse comma-separated strings to numeric vectors
  if (is.character(exon_starts)) {
    starts <- as.numeric(strsplit(as.character(exon_starts), ",")[[1]])
  } else if (is.numeric(exon_starts)) {
    starts <- as.numeric(exon_starts)
  } else {
    stop("exon_starts must be character (comma-separated) or numeric vector")
  }

  if (is.character(exon_ends)) {
    ends <- as.numeric(strsplit(as.character(exon_ends), ",")[[1]])
  } else if (is.numeric(exon_ends)) {
    ends <- as.numeric(exon_ends)
  } else {
    stop("exon_ends must be character (comma-separated) or numeric vector")
  }

  if (!strand %in% c("+", "-")) {
    stop(paste0("Invalid strand: '", strand, "'. Must be '+' or '-'"))
  }

  # Validation
  if (validate) {
    if (any(is.na(starts)) || any(is.na(ends))) {
      stop("Exon coordinates contain NA values")
    }
    if (length(starts) != length(ends)) {
      stop(paste0(
        "Mismatch: ",
        length(starts),
        " exon starts but ",
        length(ends),
        " exon ends"
      ))
    }
    if (length(starts) == 0) {
      stop("No exons provided")
    }
    if (any(starts >= ends)) {
      bad_idx <- which(starts >= ends)
      stop(paste0(
        "Invalid exon(s): start >= end at index ",
        paste(bad_idx, collapse = ", ")
      ))
    }
    if (length(starts) > 1) {
      for (i in 1:(length(starts) - 1)) {
        for (j in (i + 1):length(starts)) {
          if (!(ends[i] < starts[j] || ends[j] < starts[i])) {
            stop(paste0(
              "Exons overlap: [",
              starts[i],
              "-",
              ends[i],
              "] and [",
              starts[j],
              "-",
              ends[j],
              "]"
            ))
          }
        }
      }
    }
    if (is.unsorted(starts)) {
      warning(
        "Exon starts are not in genomic order (left->right). Reordering..."
      )
      order_idx <- order(starts)
      starts <- starts[order_idx]
      ends <- ends[order_idx]
    }
  }

  # Build exon structure in genomic order
  n_exons <- length(starts)
  exon_lengths <- ends - starts + 1

  df_genomic <- data.frame(
    genomic_idx = 1:n_exons,
    start = starts,
    end = ends,
    length = exon_lengths,
    stringsAsFactors = FALSE
  )

  # Reorder to transcript order (5'->3') based on strand
  if (strand == "+") {
    df_transcript <- df_genomic
    df_transcript$exon_num <- 1:n_exons
  } else {
    df_transcript <- df_genomic[n_exons:1, ]
    df_transcript$exon_num <- 1:n_exons
    rownames(df_transcript) <- NULL
  }

  df_transcript$is_first <- FALSE
  df_transcript$is_last <- FALSE
  df_transcript$is_first[1] <- TRUE
  df_transcript$is_last[n_exons] <- TRUE

  df_transcript <- df_transcript[, c(
    "exon_num",
    "genomic_idx",
    "start",
    "end",
    "length",
    "is_first",
    "is_last"
  )]

  structure <- list(
    exons = df_transcript,
    strand = strand,
    gene_start = min(starts),
    gene_end = max(ends),
    total_exonic_length = sum(exon_lengths),
    n_exons = n_exons,
    is_single_exon = (n_exons == 1)
  )
  class(structure) <- c("exon_structure", "list")
  return(structure)
}

# ---- get_first_exon --------------------------------------------------------
get_first_exon <- function(exon_structure) {
  if (is.null(exon_structure)) {
    return(NULL)
  }
  first <- exon_structure$exons[1, ]
  list(
    start = first$start,
    end = first$end,
    length = first$length,
    exon_num = 1
  )
}

# ---- get_last_exon ---------------------------------------------------------
get_last_exon <- function(exon_structure) {
  if (is.null(exon_structure)) {
    return(NULL)
  }
  n <- exon_structure$n_exons
  last <- exon_structure$exons[n, ]
  list(start = last$start, end = last$end, length = last$length, exon_num = n)
}

# ---- position_in_exon ------------------------------------------------------
# Check if a genomic position falls within any exon. Returns list: in_exon (logical), exon_num (transcript order), genomic_idx.
position_in_exon <- function(position, exon_structure) {
  if (is.null(exon_structure) || is.na(position)) {
    return(list(
      in_exon = FALSE,
      exon_num = NA_integer_,
      genomic_idx = NA_integer_
    ))
  }
  exons <- exon_structure$exons
  for (i in 1:nrow(exons)) {
    if (position >= exons$start[i] && position <= exons$end[i]) {
      return(list(
        in_exon = TRUE,
        exon_num = exons$exon_num[i],
        genomic_idx = exons$genomic_idx[i]
      ))
    }
  }
  return(list(
    in_exon = FALSE,
    exon_num = NA_integer_,
    genomic_idx = NA_integer_
  ))
}

# ---- position_in_first_exon ------------------------------------------------
position_in_first_exon <- function(position, exon_structure) {
  result <- position_in_exon(position, exon_structure)
  return(result$in_exon && result$exon_num == 1)
}

# ---- position_is_intronic --------------------------------------------------
# A position is intronic if it is within gene boundaries but NOT in any exon.
position_is_intronic <- function(position, exon_structure) {
  if (is.null(exon_structure) || is.na(position)) {
    return(FALSE)
  }
  within_gene <- position >= exon_structure$gene_start &&
    position <= exon_structure$gene_end
  if (!within_gene) {
    return(FALSE)
  }
  in_exon_result <- position_in_exon(position, exon_structure)
  return(!in_exon_result$in_exon)
}

# ---- get_first_intron ------------------------------------------------------
# Returns boundaries of the first intron (between first and second exon in transcript order). Returns NULL for single-exon genes.
get_first_intron <- function(exon_structure) {
  if (is.null(exon_structure) || exon_structure$is_single_exon) {
    return(NULL)
  }
  exons <- exon_structure$exons
  strand <- exon_structure$strand
  first_exon <- exons[exons$exon_num == 1, ]
  second_exon <- exons[exons$exon_num == 2, ]
  if (strand == "+") {
    intron_start <- first_exon$end + 1
    intron_end <- second_exon$start - 1
  } else {
    intron_start <- second_exon$end + 1
    intron_end <- first_exon$start - 1
  }
  list(
    start = intron_start,
    end = intron_end,
    length = intron_end - intron_start + 1
  )
}

# ---- position_in_first_intron ----------------------------------------------
position_in_first_intron <- function(position, exon_structure) {
  intron <- get_first_intron(exon_structure)
  if (is.null(intron) || is.na(position)) {
    return(FALSE)
  }
  return(position >= intron$start && position <= intron$end)
}

# ---- find_intron_number ----------------------------------------------------
# Find which intron (transcript order) contains a given position. Returns NA if position is not intronic or gene is single-exon.
find_intron_number <- function(position, exon_structure) {
  if (is.null(exon_structure) || is.na(position)) {
    return(NA_integer_)
  }
  if (exon_structure$is_single_exon) {
    return(NA_integer_)
  }
  if (!position_is_intronic(position, exon_structure)) {
    return(NA_integer_)
  }

  exons <- exon_structure$exons
  n_exons <- exon_structure$n_exons
  strand <- exon_structure$strand

  for (i in 1:(n_exons - 1)) {
    exon_i <- exons[exons$genomic_idx == i, ]
    exon_next <- exons[exons$genomic_idx == i + 1, ]
    if (position > exon_i$end && position < exon_next$start) {
      if (strand == "+") {
        return(as.integer(i))
      } else {
        return(as.integer(n_exons - i))
      }
    }
  }
  warning(sprintf(
    "Position %d is intronic but no intron found. This is a bug.",
    position
  ))
  return(NA_integer_)
}

# ---- calc_cumulative_exonic_distance ---------------------------------------
# Calculate cumulative exonic base pairs from TSS to a position, counting exonic bases only (introns skipped). Returns NA if position is not in any exon.
calc_cumulative_exonic_distance <- function(position, exon_structure) {
  if (is.null(exon_structure) || is.na(position)) {
    return(NA_real_)
  }
  pos_info <- position_in_exon(position, exon_structure)
  if (!pos_info$in_exon) {
    return(NA_real_)
  }

  exons <- exon_structure$exons
  strand <- exon_structure$strand
  position_exon_num <- pos_info$exon_num

  cumulative <- 0
  if (position_exon_num > 1) {
    for (i in 1:(position_exon_num - 1)) {
      cumulative <- cumulative + exons$length[i]
    }
  }

  position_exon <- exons[exons$exon_num == position_exon_num, ]
  if (strand == "+") {
    partial <- position - position_exon$start + 1
  } else {
    partial <- position_exon$end - position + 1
  }
  cumulative <- cumulative + partial
  return(as.numeric(cumulative))
}

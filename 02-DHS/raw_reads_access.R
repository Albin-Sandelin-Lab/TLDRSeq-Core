# raw_reads_access.R -- load TLDR-seq BED12 reads for DHS regions and samples

# ---- safe_read_table -------------------------------------------------------
# Wrap read.table() so that an empty file (no reads for this DHS x sample) returns an empty data.frame instead of raising an error. Non-empty errors are re-raised, so a real problem is never hidden.
#
# Parameters:
#   file — character, path to the table
#   ...  — passed through to read.table()
#
# Returns:
#   data.frame; zero rows if the file has no lines.
safe_read_table <- function(file, ...) {
  tryCatch(
    {
      data <- read.table(file, ...)
      return(data)
    },
    error = function(e) {
      if (grepl("no lines available in input", e$message)) {
        return(data.frame())
      } else {
        stop(e)
      }
    }
  )
}

# ---- read_TLDR_reads_of_one_DHS_in_one_sample ------------------------------
# Load the BED12 reads file for a single DHS in a single replicate.
#
# File layout on disk, one directory per DHS x sample:
#   <dir_files>/<DHS>/<sample>/reads_score_by_pA_length.bed
#   <dir_files>/<DHS>/<sample>/empty          (marker file: DHS not expressed)
#   <dir_files>/<DHS>/<sample>/done           (marker file: quantification done)
#
# Parameters:
#   query_DHS       — character, single DHS id e.g. "chr1_12345-12500"
#   query_replicate — character, single sample id e.g. "mES_WT_rep1"
#   dir_files       — directory holding one subdirectory per DHS
#
# Returns:
#   data.frame with the columns V1-V14 in BED12+extra format, or an empty data.frame if the DHS is not expressed in that sample.
#
# Column meaning (BED12 extended):
#   V1  = chromosome           V8  = thickEnd (3' end of aligned region)
#   V2  = start (0-based)      V9  = colorRGB
#   V3  = end                  V10 = number of exons in this read
#   V4  = read ID              V11 = exon sizes (comma-separated)
#   V5  = polyA tail length    V12 = exon starts (comma-separated, relative to V2)
#   V6  = strand (+/-)         V13 = DHS id
#   V7  = thickStart (5' end)  V14 = sample id
read_TLDR_reads_of_one_DHS_in_one_sample <- function(
  query_DHS,
  query_replicate,
  dir_files
) {
  stopifnot(length(query_DHS) == 1)
  stopifnot(length(query_replicate) == 1)

  tmp_dir_dhs <- file.path(dir_files, query_DHS, query_replicate)
  tmp_file_reads <- file.path(tmp_dir_dhs, "reads_score_by_pA_length.bed")

  if (!file.exists(tmp_file_reads)) {
    # "empty" marker = DHS not expressed in this sample (expected, not an error)
    if (file.exists(file.path(tmp_dir_dhs, "empty"))) {
      return(data.frame())
    } else {
      stop(sprintf(
        "DHS %s is not quantified in sample %s.\nNeither BED file nor 'empty' marker found in: %s",
        query_DHS,
        query_replicate,
        tmp_dir_dhs
      ))
    }
  } else {
    return(safe_read_table(
      tmp_file_reads,
      colClasses = c(
        "character",
        "integer",
        "integer",
        "character",
        "integer",
        "character",
        "integer",
        "integer",
        "character",
        "integer",
        "character",
        "character",
        "character",
        "character"
      )
    ))
  }
}

# ---- read_TLDR_reads_by_list_of_DHS ----------------------------------------
# Batch-load reads for all combinations of a DHS list x replicate list. Uses pbapply::pblapply so long batches show progress.
#
# Parameters:
#   query_list_DHS         — character vector of DHS ids
#   query_list_replicates  — character vector of sample ids
#   dir_files              — directory holding one subdirectory per DHS
#
# Returns:
#   data.frame with columns V1-V14, all reads stacked. Empty DHS x sample combinations contribute zero rows.
read_TLDR_reads_by_list_of_DHS <- function(
  query_list_DHS,
  query_list_replicates,
  dir_files
) {
  combinations <- expand.grid(
    DHS = unique(query_list_DHS),
    replicate = unique(query_list_replicates)
  )
  arg_list <- split(combinations, seq(nrow(combinations)))

  tmp_list <- pblapply(arg_list, function(args) {
    read_TLDR_reads_of_one_DHS_in_one_sample(
      args$DHS,
      args$replicate,
      dir_files = dir_files
    )
  })

  return(bind_rows(tmp_list))
}

# ---- rename_reads_columns --------------------------------------------------
# Rename the opaque V1-V14 columns to human-readable names. Call this after read_TLDR_reads_by_list_of_DHS().
#
# Parameters:
#   df_reads — data.frame with columns V1-V14 (output of read_TLDR_reads_*)
#
# Returns:
#   data.frame with named columns: chromosome, start, end, read_id, polyA_tail_length, strand, thickStart, thickEnd, colorRGB, num_of_exons, exon_sizes, exon_starts, DHS, sample
rename_reads_columns <- function(df_reads) {
  df_reads %>%
    dplyr::rename(
      chromosome = V1,
      start = V2,
      end = V3,
      read_id = V4,
      polyA_tail_length = V5,
      strand = V6,
      thickStart = V7,
      thickEnd = V8,
      colorRGB = V9,
      num_of_exons = V10,
      exon_sizes = V11,
      exon_starts = V12,
      DHS = V13,
      sample = V14
    )
}

# ---- read_total_depth_of_one_DHS -------------------------------------------
# Read the total read count for a single DHS x sample from the "done" marker.
#
# Parameters:
#   dir_files — directory holding one subdirectory per DHS
#   DHS       — character, single DHS id
#   sampleID  — character, single sample id
#
# Returns:
#   numeric total count; 0 if the "empty" marker exists; NA if neither marker is found.
read_total_depth_of_one_DHS <- function(dir_files, DHS, sampleID) {
  tmp_dir_dhs <- file.path(dir_files, DHS, sampleID)
  if (!file.exists(file.path(tmp_dir_dhs, "done"))) {
    if (file.exists(file.path(tmp_dir_dhs, "empty"))) {
      return(0)
    } else {
      return(NA)
    }
  } else {
    return(read.table(file.path(tmp_dir_dhs, "done"))$V1)
  }
}

# ---- batch_summarise_total_depth_of_DHS ------------------------------------
# Summarise total read counts across a DHS list x sample list.
#
# Parameters:
#   dir_files   — directory holding one subdirectory per DHS
#   DHS_list    — character vector of DHS ids
#   sample_list — character vector of sample ids
#
# Returns:
#   data.frame with columns DHS, sample, total_count.
batch_summarise_total_depth_of_DHS <- function(
  dir_files,
  DHS_list,
  sample_list
) {
  combinations <- expand.grid(
    DHS = unique(DHS_list),
    replicate = unique(sample_list)
  )
  arg_list <- split(combinations, seq(nrow(combinations)))

  tmp_list <- lapply(arg_list, function(args) {
    read_total_depth_of_one_DHS(
      dir_files = dir_files,
      DHS = args$DHS,
      sampleID = args$replicate
    )
  })

  return(data.frame(
    DHS = combinations$DHS,
    sample = combinations$replicate,
    total_count = unlist(tmp_list)
  ))
}

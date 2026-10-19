# config_and_libraries.R -- shared libraries, sample universe, input-path convention

# ---- Libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(GenomicRanges)
  library(S4Vectors)
  library(pbapply)
})

# ---- Sample universe ----
# The samples cover wild type and two knockdowns: WT, RBM7 and ZCCHC8.
all_replicates <- c(
  "mES_WT_rep1",
  "mES_WT_rep2",
  "mES_WT_rep3",
  "mES_WT_rep4",
  "mES_WT_rep5",
  "mES_WT_rep6",
  "mES_RBM7_rep1",
  "mES_RBM7_rep2",
  "mES_RBM7_rep3",
  "mES_ZCCHC8_rep1",
  "mES_ZCCHC8_rep2",
  "mES_ZCCHC8_rep3"
)

# ---- env_path ----
#
# Parameters:
#   name — character, name of the environment variable that holds the path
#
# Returns:
#   character, the path. The call stops, naming the variable, when it is unset or empty.
env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop(
      sprintf("Environment variable %s must name the file.", name),
      call. = FALSE
    )
  }
  value
}

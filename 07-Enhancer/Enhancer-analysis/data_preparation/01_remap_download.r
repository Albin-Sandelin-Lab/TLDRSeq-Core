## Download the ReMap 2022 mouse ESC bigWig archive (mm10) from Recherche Data Gouv
## (doi:10.57745/GEIHVN): 7 parts of one tar archive.
##
## Usage (from 07-Enhancer/Enhancer-analysis/data_preparation/): Rscript 01_remap_download.r
## Output: $REMAP_DIR/ReMap_mmus_mm10_mESC_archive_part.tar.part.a{a..g}
##         (the metadata file ReMap_Mmus_mESC_metafile.txt is part of the same dataset)

library(parallel)
library(dataverse)

source("../../../config/paths.R")

Sys.setenv("DATAVERSE_SERVER" = "entrepot.recherche.data.gouv.fr")

dir.create(REMAP_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(REMAP_DIR)

basename <- "ReMap_mmus_mm10_mESC_archive_part.tar.part.a"
index <- c("a", "b", "c", "d", "e", "f", "g")
files.names <- paste0(basename, index)

download_results <- mclapply(files.names, function(X) {
  # one failed download doesn't stop the others
  tryCatch({
    raw_bytes <- get_file_by_name(
      filename = X,
      dataset = "doi:10.57745/GEIHVN"
    )
    writeBin(raw_bytes, X)
    return(paste("Successfully downloaded and saved:", X))
  }, error = function(e) {
    return(paste("Error downloading", X, ":", e$message))
  })
}, mc.cores = 10)

cat(paste(unlist(download_results), collapse = "\n"))

## Select the ReMap mESC experiments in control / wild-type conditions and pool them per TF.
##
## Usage (from 07-Enhancer/Enhancer-analysis/data_preparation/): Rscript 03_remap_select_pool.r
##
## Kept: experiments without condition or with a control / baseline condition (WT, parental,
## control knock-down, vehicle, 0 h / d0, untreated, 2i, serum); 48 h / 72 h time points removed.
## Pooled signal of a TF = sum of the coverage of its experiments.
##
## Input : $REMAP_DIR/ReMap_Mmus_mESC_metafile.txt, $REMAP_DIR/bigwig_files_mm39/ (02)
## Output: $REMAP_DIR/cleaned_mapping_table.txt   selected experiments
##         $REMAP_DIR/pooled_WT_bigWig/Pooled_<TF>.bw

library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(tidyverse)

source("../../../config/paths.R")

setwd(REMAP_DIR)


#### Control / wild-type experiments

# fill = TRUE handles rows with only 3 columns (where the condition is absent)
remap_meta <- read.delim("ReMap_Mmus_mESC_metafile.txt",
                         header = FALSE,
                         col.names = c("source_accession", "target", "biosample", "condition"),
                         fill = TRUE,
                         stringsAsFactors = FALSE)

# control / baseline conditions, case-insensitive exact matches (not "2i-d3" or "shGFP-KO")
ctrl_pattern <- "^(?i)(wt.*|parental|es-wt|sictrl|shscr|shluc|shgfp|shctr|control.*|scr|ctrl|dmso|etoh|mock.*|0h|d0|untreated|unt|asynchronous|2i|e14-2i|serum|e14-serum)$"

clean_meta <- remap_meta %>%
    mutate(condition = ifelse(is.na(condition), "", trimws(condition))) %>%
    filter(condition == "" | grepl(ctrl_pattern, condition, perl = TRUE)) %>%
    filter(!grepl("48|72", condition))

print(dim(clean_meta))
write.table(clean_meta, "cleaned_mapping_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#### Pool the experiments of each TF

meta_data <- read.table("cleaned_mapping_table.txt", sep = "\t", header = TRUE) %>%
    mutate(basename = paste0(source_accession, ".", target, ".", biosample, "_", condition, ".bigwig.bw")) %>%
    mutate(basename = basename %>% str_replace("_\\.", "."))

out_dir <- "pooled_WT_bigWig"
dir.create(out_dir, showWarnings = FALSE)

pool_bw_files <- function(current_tf, meta_data, out_dir) {
    tf_datasets <- meta_data %>%
        filter(target == current_tf) %>%
        pull(basename)
    file_names <- file.path("bigwig_files_mm39", tf_datasets)
    valid_files <- file_names[file.exists(file_names)]
    output_name <- paste0("Pooled_", current_tf, ".bw")
    if (length(valid_files) > 1) {
        ## sum of the coverage of all experiments
        gr_list <- lapply(valid_files, import.bw)
        combined_gr <- do.call(c, gr_list)
        pooled_signal <- coverage(combined_gr, weight = "score")
        export.bw(pooled_signal, file.path(out_dir, output_name))
        rm(gr_list, combined_gr, pooled_signal)
        gc()
        return(paste("Processed:", current_tf))
    } else if (length(valid_files) == 1) {
        file.copy(valid_files, file.path(out_dir, output_name), overwrite = TRUE)
        return(paste("Copied:", current_tf, "(only one file)"))
    } else {
        return(paste("Skipped:", current_tf, "(zero files found)"))
    }
}

tf_list <- unique(meta_data$target)

results <- mclapply(tf_list,
                    FUN = pool_bw_files,
                    meta_data = meta_data,
                    out_dir = out_dir,
                    mc.cores = 10)

print(unlist(results))

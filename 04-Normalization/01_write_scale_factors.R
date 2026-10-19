## chrM scale factor of each TLDRSeq_sup replicate: chrM depth / mean chrM depth over the
## 21 replicates (define_scale_factor(), scale_factors.R).
##
## Usage (from 04-Normalization/): Rscript 01_write_scale_factors.R
##
## Input : CHRM_DEPTH (config/paths.R): "<depth> <path/to/bed12_file>" per replicate
## Output: SCALE_FACTORS_CSV (config/paths.R): sample ("scale_<replicate>"), scale_factor

library(tidyverse)

source("../config/paths.R")
source("scale_factors.R")

all_replicates <- c("mES_WT_rep1", "mES_WT_rep2", "mES_WT_rep3", "mES_WT_rep4", "mES_WT_rep5", "mES_WT_rep6",
                    "mES_ARS2_rep1", "mES_ARS2_rep2", "mES_ARS2_rep3",
                    "mES_ZC3H18_rep1", "mES_ZC3H18_rep2", "mES_ZC3H18_rep3",
                    "mES_ZC3H4_rep1", "mES_ZC3H4_rep2", "mES_ZC3H4_rep3",
                    "mES_RBM7_rep1", "mES_RBM7_rep2", "mES_RBM7_rep3",
                    "mES_ZCCHC8_rep1", "mES_ZCCHC8_rep2", "mES_ZCCHC8_rep3")

df_scale_factor <- define_scale_factor(CHRM_DEPTH, sample_levels = all_replicates)
validate_scale_factor_table(df_scale_factor)

dir.create(dirname(SCALE_FACTORS_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(df_scale_factor %>% mutate(sample = paste0("scale_", sample)), SCALE_FACTORS_CSV)

## Effect of AT content on PROMPT sensitivity to RBM7 and ZCCHC8 depletion.
##
## Usage (from 05-PROMPT/): Rscript 07_prompt_AT_model.r
##
## For each PROMPT window (06_prompt_windows.r), the log2FC (RBM7 or ZCCHC8 vs WT) is modelled
## as a function of the cumulative AT count from the PROMPT TSS, correcting for the distance
## from the TSS:
##   log2FC ~ log10(cumulative AT count) + log10(cumulative width)
## The length effect is then removed (log2FC - beta_length * log10(width)) and the corrected
## log2FC is plotted against the cumulative AT count.
## Windows: median read count > 3 in the genotype; log2FC < 0.1 set to 0.
##
## Input : $PROMPT_DIR/prompt_windows/at_log2FC_per_window_per_DHS.tsv (06_prompt_windows.r)
## Output: $PROMPT_DIR/AT_model/AT_effect_prompts_{RBM7,ZCCHC8}.pdf

library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggpmisc)

source("../config/paths.R")

out_dir <- file.path(PROMPT_DIR, "AT_model")
dir.create(out_dir, showWarnings=FALSE, recursive =TRUE)

at.log2FC.prompts.df <- fread(file.path(PROMPT_DIR, "prompt_windows", "at_log2FC_per_window_per_DHS.tsv"))

custom_blues <- RColorBrewer::brewer.pal(9, "Blues")[3:9]

for (geno.elt in c("RBM7", "ZCCHC8")) {
    window.infos.clean <- at.log2FC.prompts.df %>%
        mutate(log2FC = ifelse(log2FoldChange < 0.1, 0, log2FoldChange)) %>%
        dplyr::filter(geno == geno.elt) %>%
        dplyr::filter(median_count > 3) %>%
        dplyr::filter(!is.na(log2FC),
                      !is.na(cumulative_AT.counts))

    ## 1. Fit the model
    ## log2FC = Intercept + Beta_AT * AT + Beta_Len * Len + Error
    current_model <- lm(log2FC ~ log10(cumulative_AT.counts) + log10(width.cumulative),
                        data = window.infos.clean)
    at_slope <- coef(current_model)["log10(cumulative_AT.counts)"]
    at_pval <- summary(current_model)$coefficients["log10(cumulative_AT.counts)", "Pr(>|t|)"]
    pval_label <- formatC(at_pval, format = "e", digits = 2)

    ## 2. Remove the length effect ("regress it out")
    ## Corrected log2FC = raw log2FC - length contribution
    ## This leaves behind the intercept + AT effect + residuals
    beta_length <- coef(current_model)["log10(width.cumulative)"]
    length_effect <- beta_length * log10(window.infos.clean$width.cumulative)
    window.infos.clean$log2FC_corrected <- window.infos.clean$log2FC - length_effect

    ## 3. Plot the corrected log2FC vs AT
    p <- ggplot(window.infos.clean, aes(x = log10(cumulative_AT.counts), y = log2FC_corrected)) +
        geom_hex() +
        scale_fill_gradientn(colors = custom_blues, trans = "log10") +
        ## Trend line (simple linear fit on the corrected data)
        geom_smooth(method = "lm", color = "blue", size = 1.2) +
        labs(title = paste("Source: prompts,", geno.elt),
             subtitle = paste0("AT Slope: ", signif(at_slope, 3), "  (p = ", pval_label, ")"),
             y = "Log2FC (Length Corrected)",
             x = "Cumulative AT Counts (Log10)") +
        stat_poly_eq(aes(label = after_stat(eq.label)), formula = y ~ x, label.x = "left", label.y = "top") +
        theme_minimal() +
        xlim(c(0.5,3))
    ggsave(file.path(out_dir, paste0("AT_effect_prompts_", geno.elt, ".pdf")), plot = p, width = 6, height = 4)
}

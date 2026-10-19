# cluster_genes_by_proportion.R -- cluster genes by 3' category proportion profile

# ---- cluster_genes_by_proportion -------------------------------------------
# Cluster genes on their 3' category proportion profile and label each DHS.
#
# The row ordering and the cluster labels both come from the drawn heatmap: draw() fixes the dendrogram, row_order() reads the displayed order, and cutree() at heatmap_cluster_k turns the tree into cluster numbers.
#
# Parameters:
#   df_freq_long       — data.frame with columns DHS, feature_id, geno, feature_count, geno_count, scaled_freq (from 06-Genic/proportions' calculate_pooled_counts_then_proportion()).
#   gene_map           — data.frame with columns DHS and gene_name, one row per DHS. Must already be unique on DHS.
#   heatmap_cluster_k  — integer. Number of clusters to cut the row dendrogram into (default 10; configured a priori, not derived from the data).
#
# Returns:
#   list with:
#     df_gene_groups             — data.frame with columns DHS, gene_name, gene_group ("cluster_<k>")
#     df_prop_wide_simp_clusters — data.frame, one row per DHS, with the proportion columns plus heatmap_cluster and heatmap_row_order
#     ht_prop                    — the annotated ComplexHeatmap object
cluster_genes_by_proportion <- function(
  df_freq_long,
  gene_map,
  heatmap_cluster_k = 10
) {
  df_prop_wide <- df_freq_long %>%
    mutate(geno_feature = paste0(geno, "_", feature_id)) %>%
    pivot_wider(
      id_cols = "DHS",
      names_from = "geno_feature",
      values_from = "scaled_freq",
      values_fill = 0
    ) %>%
    column_to_rownames("DHS")
  col_fun_prop <- circlize::colorRamp2(
    c(0, 1), # numeric breakpoints
    c("white", "#028090") # corresponding colors
  )
  # keep the proportion columns of the three genotyped samples
  df_prop_wide_simp <- df_prop_wide %>%
    dplyr::select(contains("WT"), contains("RBM7"), contains("ZCCHC8"))

  ht_prop <- df_prop_wide_simp %>%
    ComplexHeatmap::pheatmap(
      show_rownames = F,
      use_raster = FALSE,
      border_color = NA,
      scale = "none",
      cluster_rows = T,
      color = col_fun_prop,
      annotation_names_row = F,
      heatmap_legend_param = list(title = "Proportion")
    )

  ht_prop_drawn <- draw(ht_prop)
  heatmap_row_order <- unlist(row_order(ht_prop_drawn), use.names = FALSE)
  heatmap_clusters <- cutree(
    as.hclust(row_dend(ht_prop_drawn)),
    k = heatmap_cluster_k
  )

  df_prop_wide_simp_clusters <- df_prop_wide_simp %>%
    as.data.frame() %>%
    tibble::rownames_to_column("DHS") %>%
    mutate(
      heatmap_cluster = unname(heatmap_clusters[DHS]),
      heatmap_row_order = match(
        DHS,
        rownames(df_prop_wide_simp)[heatmap_row_order]
      )
    ) %>%
    arrange(heatmap_row_order)

  gene_map <- gene_map %>%
    distinct(DHS, gene_name)

  if (nrow(gene_map) != n_distinct(gene_map$DHS)) {
    stop("gene_map must contain exactly one gene_name per DHS.")
  }

  df_gene_groups <- df_prop_wide_simp_clusters %>%
    transmute(
      DHS,
      gene_group = paste0("cluster_", heatmap_cluster)
    ) %>%
    left_join(gene_map, by = "DHS") %>%
    select(DHS, gene_name, gene_group)

  if (
    any(is.na(df_gene_groups$gene_name) | !nzchar(df_gene_groups$gene_name))
  ) {
    stop("Missing gene_name after joining heatmap clusters to the gene map.")
  }

  ht_prop <- df_prop_wide_simp %>%
    ComplexHeatmap::pheatmap(
      show_rownames = F,
      border_color = NA,
      scale = "none",
      cluster_rows = T,
      color = col_fun_prop,
      annotation_names_row = F,
      annotation_row = df_prop_wide_simp_clusters %>%
        transmute(DHS, cluster = factor(heatmap_cluster)) %>%
        column_to_rownames("DHS"),
      annotation_colors = list(
        cluster = c(
          "1" = "#1b9e77",
          "2" = "#d95f02",
          "3" = "#7570b3",
          "4" = "#e7298a",
          "5" = "#66a61e",
          "6" = "#e6ab02",
          "7" = "#a6761d",
          "8" = "#666666",
          "9" = "#1f78b4",
          "10" = "#33a02c"
        )
      ),
      heatmap_legend_param = list(title = "Proportion")
    )

  return(list(
    df_gene_groups = df_gene_groups,
    df_prop_wide_simp_clusters = df_prop_wide_simp_clusters,
    ht_prop = ht_prop
  ))
}

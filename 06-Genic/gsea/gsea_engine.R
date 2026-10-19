# gsea_engine.R -- the GSEA engine functions gsea_of_proportions.Rmd calls
#
# The implementations follow the myExpr package, reduced to the functions the notebook reaches.

# ---- complete_id_map: used by gsea_by_one_coef
complete_id_map = function(df, gene_id_map){

  if("entrez_id" %in% colnames(df)){
    if("gene_id" %in% colnames(df)){
      out = df
    }else{
      out = df %>%
        left_join(as.data.frame(gene_id_map), by = "entrez_id")
    }

  }else if("gene_id" %in% colnames(df)){
    out = df %>%
      left_join(as.data.frame(gene_id_map), by = "gene_id")
  }else{
    stop("gene_id or entrez_id not in colname(df)")
  }

  return(out)
}

# ---- gsea_by_one_coef: ranks genes and runs the GO GSEA
gsea_by_one_coef = function(df_coef, coef_name, gene_id_map,
                            pcutoff = 0.05,
                            GO_or_KEGG = "GO",
                            species = "mouse",
                            custom_TERM2GENE = NA,
                            custom_TERM2NAME = NA,
                            ...){

  df_coef = complete_id_map(df_coef, gene_id_map) %>% 
    as.data.frame()

  if(species == "mouse"){
    library(org.Mm.eg.db)
    orgdb = org.Mm.eg.db
    species_wp = "Mus musculus"
  }else if(species == "human"){
    library(org.Hs.eg.db)
    orgdb = org.Hs.eg.db
    species_wp = "Homo sapiens"
  }else{
    stop("Species Not Supported")
  }

  df_coef_clean = df_coef %>%
    filter(!duplicated(entrez_id))

  geneList = df_coef_clean[,coef_name]
  
  names(geneList) = df_coef_clean$entrez_id
  geneList = sort(geneList, decreasing = TRUE)

  if(GO_or_KEGG == "GO"){
    outout = clusterProfiler::gseGO(geneList,
                                    OrgDb = orgdb,ont = "BP",

                                    pvalueCutoff = pcutoff,
                                    ...)
  }else if(GO_or_KEGG == "KEGG"){
    outout = clusterProfiler::gseKEGG(geneList,
                                      organism = species,
                                      pvalueCutoff = pcutoff,
                                      ...)
  }else if(GO_or_KEGG == "MKEGG"){
    outout = clusterProfiler::gseMKEGG(geneList,
                                      organism = species,
                                      pvalueCutoff = pcutoff,
                                      ...)
    
  }else if(GO_or_KEGG == "GO_ALL"){
    outout = clusterProfiler::gseGO(geneList,
                                    OrgDb = orgdb,ont = "ALL",

                                    pvalueCutoff = pcutoff,
                                    ...)
  }else if(GO_or_KEGG == "custom" &
           length(custom_TERM2GENE) > 1 ){
    outout = clusterProfiler::GSEA(geneList,
                                   minGSSize = 1,maxGSSize = Inf,
                                   TERM2GENE = custom_TERM2GENE,
                                   TERM2NAME = custom_TERM2NAME,
                                   nPermSimple = 10000,
                                    pvalueCutoff = pcutoff)
  }else if(GO_or_KEGG == "WikiPathway"){
    outout = clusterProfiler::gseWP(geneList,
                                    organism = species_wp,
                                    pvalueCutoff = pcutoff)
    
  }else{

    print(paste0("Unsupported GO_or_KEGG: ", GO_or_KEGG))
    }

  return(outout)
}

# ---- fortify.compareClusterResult.withfilter: produces the tileplot data frame
fortify.compareClusterResult.withfilter <- function(model, data, showCategory=5,
                                         by="geneRatio",
                                         split=NULL,
                                         includeAll=TRUE,
                                         excludeCategory = c(),
                                         showCertainCategory = c()) {
  clProf.df <- as.data.frame(model)
  # exclude everything in excludeCategory, a vector of terms
  if(length(showCertainCategory) == 0){
    clProf.df = clProf.df[!(clProf.df$Description %in% excludeCategory),]
  }else{
    clProf.df = clProf.df[(clProf.df$Description %in% showCertainCategory),]
  }

  .split <- split

  ## get top 5 (default) categories of each gene cluster.
  if (is.null(showCategory)) {
    result <- clProf.df
  } else {
    Cluster <- NULL # to satisfy codetools

    topN <- function(res, showCategory) {
      plyr::ddply(.data = res,
            .variables = plyr::`.`(Cluster),
            .fun = function(df, N) {
              if (length(df$Count) > N) {
                if (any(colnames(df) == "pvalue")) {
                  idx <- order(df$pvalue, decreasing=FALSE)[1:N]
                } else {
                  ## for groupGO
                  idx <- order(df$Count, decreasing=T)[1:N]
                }
                return(df[idx,])
              } else {
                return(df)
              }
            },
            N=showCategory
      )

    }

    if (!is.null(.split) && .split %in% colnames(clProf.df)) {
      lres <- split(clProf.df, as.character(clProf.df[, .split]))
      lres <- lapply(lres, topN, showCategory = showCategory)
      result <- do.call('rbind', lres)
    } else {
      result <- topN(clProf.df, showCategory)
    }

  }

  ID <- NULL
  if (includeAll == TRUE) {
    result <- subset(clProf.df, ID %in% result$ID)
  }

  ## remove zero count
  result$Description <- as.character(result$Description) ## un-factor
  GOlevel <- result[,c("ID", "Description")] ## GO ID and Term
  GOlevel <- unique(GOlevel)

  result <- result[result$Count != 0, ]
  result$Description <- factor(result$Description,
                               levels=rev(GOlevel[,2]))


  if (by=="rowPercentage") {
    Description <- Count <- NULL # to satisfy codetools
    result <- ddply(result,
                    .(Description),
                    transform,
                    Percentage = Count/sum(Count),
                    Total = sum(Count))

    ## label GO Description with gene counts.
    x <- mdply(result[, c("Description", "Total")], paste, sep=" (")
    y <- sapply(x[,3], paste, ")", sep="")
    result$Description <- y

    ## restore the original order of GO Description
    xx <- result[,c(2,3)]
    xx <- unique(xx)
    rownames(xx) <- xx[,1]
    Termlevel <- xx[as.character(GOlevel[,1]),2]

    ##drop the *Total* column
    result <- result[, colnames(result) != "Total"]

    result$Description <- factor(result$Description,
                                 levels=rev(Termlevel))

  } else if (by == "count") {
    ## nothing
  } else if (by == "geneRatio") {
    gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
    gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
    result$GeneRatio <- gsize/gcsize
    cluster <- paste(as.character(result$Cluster),"\n", "(", gcsize, ")",
                     sep="")
    lv <- unique(cluster)[order(as.numeric(unique(result$Cluster)))]
    result$Cluster <- factor(cluster, levels = lv)
  } else {
    ## nothing
  }
  return(result)
}

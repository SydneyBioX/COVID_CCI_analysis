

pCCI_byCellType <- function(object) {
  object1 <- methods::slot(object, "netP")
  prob1 <- object1$prob
  df <- melt(apply(prob1, 3, function(x) {
    df <- melt(x)
    colnames(df) <- c("Ligand", "Receptor", "value")
    df
  }))
  df <- df[, c("Ligand", "Receptor", "L1", "value")]
  colnames(df)[3] <- "Pathway"
  return(df)
}

pCCI_summarise <- function(pCCI_per_sample) {
  df <- melt(pCCI_per_sample)
  res <- aggregate(df$value, 
                   list(df$Ligand,
                        df$Receptor,
                        df$L1,
                        df$Pathway),
                   sum)
  colnames(res) <- c("Ligand", 
                     "Receptor",
                     "sample",
                     "Pathway",
                     "value")
  features <- paste(res$Ligand,
                    res$Receptor,
                    res$Pathway, sep = "_")
  
  res$features <- features
  
  res$value <- (res$value  - min(res$value ))/(max(res$value ) - min(res$value ))
  
  return(res)
}


plot_pathwayCCI <- function(pCCI_all, ligand, receptor, nclust_pathway = 4, anno_col = NULL) {
  keep <- pCCI_all$Ligand %in% ligand &
    pCCI_all$Receptor %in% receptor
  pmat <- pCCI_all[keep, ] %>%
    dcast2(Pathway~sample, 
           fun.aggregate = sum, value.var = "value")
  pmat <- pmat[rowSums(pmat) != 0 & rowSums(pmat != 0) > 3, ]
  
  
  hclust_pathway <- hclust(dist(t(apply(-1/log10(pmat), 1, scale))))
  pathway_clust <- cutree(hclust_pathway, k = nclust_pathway)
  
  anno_row <- data.frame(pathway_cluster = factor(pathway_clust))
  
  rownames(anno_row) <- names(pathway_clust)
  anno_color <- list()
  anno_color$pathway_cluster <- RColorBrewer::brewer.pal(nclust_pathway, "Set2")
  names(anno_color$pathway_cluster) <- seq_len(nclust_pathway)
  pheatmap(-1/log10(pmat), 
           scale = "row",
           breaks = seq(-2, 2, len = 100),
           annotation_col = data.frame(anno_col),
           annotation_row = anno_row,
           main = paste("Source:", ligand, "- Target:", receptor),
           annotation_colors = anno_color)
}



selectCCIfeatures <- function(pCCI_all_mat, groupings) {
  kruskal_pvalue <- list()
  for (i in 1:nrow(pCCI_all_mat)) {
    #if (i %% 100 == 0) cat(i, "...")
    kruskal_res <- try(kruskal.test(unlist(pCCI_all_mat[i,]) ~ groupings), silent = TRUE)
    kruskal_pvalue[[i]] <- try(kruskal_res$p.value, silent = TRUE)
    
  }
  
  kruskal_pvalue <- lapply(kruskal_pvalue, function(x) {
    if (class(x) == "try-error") {
      x <- NULL
    }
    x
  })
  names(kruskal_pvalue) <- rownames(pCCI_all_mat)
  kruskal_pvalue <- unlist(kruskal_pvalue)
  
  kruskal_pvalue <- p.adjust(kruskal_pvalue, method = "BH")
  return(kruskal_pvalue)
}

tCCI <- function(pCCI_per_sample) {
  aff_mat_bySample <- lapply(pCCI_per_sample,
                             function(x) dcast2(x, Ligand~Receptor,
                                                fun.aggregate = sum, value.var = "value"))
  
  aff_mat_bySample <- lapply(aff_mat_bySample, function(x) {
    mat <- matrix(0, ncol = length(all_cellTypes), nrow = length(all_cellTypes))
    colnames(mat) <- rownames(mat) <- all_cellTypes
    mat[rownames(x), colnames(x)] <- as.matrix(x)
    mat
  })
  
  aff_mat_bySample <- lapply(aff_mat_bySample, function(x) {
    (x - min(x))/(max(x) - min(x))
  })
  return(aff_mat_bySample)
}


plotTCCI <- function(tCCI_results, title = NULL) {
  
  pheatmap(tCCI_results, cluster_cols = FALSE, 
           cluster_rows = FALSE,
           main = title,
           color =  colorRampPalette(c("white", 
                                       brewer.pal(n = 7, 
                                                  name = "Reds")))(100),
           breaks = seq(0, max(tCCI_results), max(tCCI_results)/100))
}


dcast2 <- function (x, ...) {
  data <- reshape2::dcast(x, ...)
  rownames(data) <- data[, 1]
  data <- data[, -1]
  return(data)
}

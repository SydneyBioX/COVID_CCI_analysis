
library(scClassify)
library(SingleCellExperiment)
library(pbmcapply)


runMultiLevelscClassify <- function(sce, trainObject, 
                                    n_cells_subset = 30000,
                                    ncores = 10) {
  
  set.seed(2020)
  #sce <- sce[rowSums(logcounts(sce)) >= 10, ]
  hvg <- unlist(trainObject$hierarchyKNNRes$limma$hvg)
  sce <- sce[rownames(sce) %in% hvg, ]
  sce
  gc(reset = TRUE)
  
  index_list <- split(seq(ncol(sce)), (seq(ncol(sce))-1) %/% n_cells_subset) 
  
  cat("Divided the number of test data into", length(index_list), "blocks \n")
  
  print("Running level 1")
  
  pred_res_level1 <- pbmclapply(1:length(index_list), function(i) {
    pred_res <- predict_scClassify(exprsMat_test = (logcounts(sce[, index_list[[i]]])),
                                   trainObject,
                                   verbose = TRUE)
    
    pred_res <- pred_res$pearson_WKNN_limma$predRes
    pred_length <- unlist(lapply(strsplit(pred_res, "_"), length))
    pred_res[pred_length >= 2] <- "intermediate"
    return(pred_res)
  }, mc.cores = ncores)
  
  pred_res_level1 <- unlist(pred_res_level1)
  pred_res_level1 <- pred_res_level1[colnames(sce)]
  
  
  print(table(pred_res_level1))
  
  
  train_idx <- !pred_res_level1 %in% c("unassigned", "intermediate", 
                                       names(table(pred_res_level1)[table(pred_res_level1) < 5]))
  
  print("Running level 2")
  
  print("Train/test split")
  print(table(train_idx))
  exprsMat_train <- logcounts(sce)[, train_idx]
  cellType_train <- pred_res_level1[train_idx]
  names(cellType_train) <- colnames(sce)[train_idx]
  exprsMat_test <- logcounts(sce)[, !train_idx]
  
  dim(exprsMat_train)
  dim(exprsMat_test)
  
  celltype_list <- names(table(cellType_train))
  
  set.seed(2020)
  train_index_list <- lapply(1:10, function(i) unlist(sapply(celltype_list, function(x) {
    idx <- which(cellType_train == x)
    if (length(idx) < 500) {
      n <- length(idx)
    } else {
      n <- max(length(idx) * 0.05, 500)
    }
    sample(idx, n)
  })))
  
  lapply(train_index_list, function(x) table(cellType_train[x]))
  
  if (ncol(exprsMat_test) > 50000) {
    test_index_list <- split(seq(ncol(exprsMat_test)), (seq(ncol(exprsMat_test))-1) %/% n_cells_subset) 
    exprsMat_test_list <- lapply(test_index_list, function(x) exprsMat_test[, x])
  } else {
    exprsMat_test_list <- list(test = exprsMat_test)
  }
  
  names(exprsMat_test_list) <- paste("test", names(exprsMat_test_list), sep = "_")
  pred_res_level2 <- pbmclapply(1:length(train_index_list), function(i) {
    trainCohort <- scClassify(exprsMat_train = as(exprsMat_train[, train_index_list[[i]]], "dgCMatrix"),
                              cellTypes_train = cellType_train[train_index_list[[i]]],
                              exprsMat_test = exprsMat_test_list,
                              similarity = c("pearson"),
                              selectFeatures = c("limma"),
                              algorithm = c("WKNN"),
                              parallel = F,
                              prob_threshold = 0.7,
                              hopach_kmax = 9,
                              verbose = T,
                              k = 5,
                              topN = 30)
    pred_res <- unlist(lapply(trainCohort$testRes, function(x) x$pearson_WKNN_limma$predRes))
    pred_length <- unlist(lapply(strsplit(pred_res, "_"), length))
    pred_res[pred_length >= 2] <- "intermediate"
    
    return(pred_res)
  }, mc.cores = ncores)
  
  pred_res_level2 <- do.call(cbind, pred_res_level2)
  
  pred_res_level2_list <- apply(pred_res_level2, 1, function(x) {
    tab <- table(x)
    c(names(tab)[which.max(tab)], tab[which.max(tab)])
  })
  pred_res_level2_list <- t(pred_res_level2_list)
  pred_res_level2_list <- data.frame(pred_res_level2_list)
  pred_res_level2_list$X2 <- as.numeric(pred_res_level2_list$X2)/10
  
  pred_res_level2_list$final <- pred_res_level2_list$X1
  pred_res_level2_list$final[pred_res_level2_list$X2 < 0.8] <- "unassigned"
  rownames(pred_res_level2_list) <- unlist(lapply(strsplit(rownames(pred_res_level2_list), "\\."), function(x) {
    x <- paste(x[-1], collapse  = ".")
    x
  }))
  
  
  pred_res_level2 <- pred_res_level1
  
  pred_res_level2[rownames(pred_res_level2_list)] <- pred_res_level2_list$final
  print(table(pred_res_level2))
  
  print("Running level 3")
  
  
  train_idx <- !pred_res_level2 %in% c("unassigned", "intermediate")
  print("Train/test split")
  print(table(train_idx))
  exprsMat_train <- logcounts(sce)[, train_idx]
  cellType_train <- pred_res_level2[train_idx]
  names(cellType_train) <- colnames(sce)[train_idx]
  exprsMat_test <- logcounts(sce)[, !train_idx]
  
  dim(exprsMat_train)
  dim(exprsMat_test)
  
  celltype_list <- names(table(cellType_train))
  
  set.seed(2020)
  train_index_list <- lapply(1:10, function(i) unlist(sapply(celltype_list, function(x) {
    idx <- which(cellType_train == x)
    if (length(idx) < 500) {
      n <- length(idx)
    } else {
      n <- max(length(idx) * 0.05, 500)
    }
    sample(idx, n)
  })))
  
  lapply(train_index_list, function(x) table(cellType_train[x]))
  
  
  
  if (ncol(exprsMat_test) > 50000) {
    test_index_list <- split(seq(ncol(exprsMat_test)), (seq(ncol(exprsMat_test))-1) %/% n_cells_subset) 
    exprsMat_test_list <- lapply(test_index_list, function(x) exprsMat_test[, x])
  } else {
    exprsMat_test_list <- list(test = exprsMat_test)
  }
  
  names(exprsMat_test_list) <- paste("test", names(exprsMat_test_list), sep = "_")
  pred_res_level3 <- pbmclapply(1:length(train_index_list), function(i) {
    trainCohort <- scClassify(exprsMat_train = as(exprsMat_train[, train_index_list[[i]]], "dgCMatrix"),
                              cellTypes_train = cellType_train[train_index_list[[i]]],
                              exprsMat_test = exprsMat_test_list,
                              similarity = c("pearson"),
                              selectFeatures = c("limma"),
                              algorithm = c("WKNN"),
                              parallel = F,
                              prob_threshold = 0.7,
                              hopach_kmax = 9,
                              verbose = T,
                              k = 5,
                              topN = 30)
    pred_res <- unlist(lapply(trainCohort$testRes, function(x) x$pearson_WKNN_limma$predRes))
    pred_length <- unlist(lapply(strsplit(pred_res, "_"), length))
    pred_res[pred_length >= 2] <- "intermediate"
    
    return(pred_res)
  }, mc.cores = ncores)
  
  pred_res_level3 <- do.call(cbind, pred_res_level3)
  
  pred_res_level3_list <- apply(pred_res_level3, 1, function(x) {
    tab <- table(x)
    c(names(tab)[which.max(tab)], tab[which.max(tab)])
  })
  pred_res_level3_list <- t(pred_res_level3_list)
  pred_res_level3_list <- data.frame(pred_res_level3_list)
  pred_res_level3_list$X2 <- as.numeric(pred_res_level3_list$X2)/10
  
  pred_res_level3_list$final <- pred_res_level3_list$X1
  pred_res_level3_list$final[pred_res_level3_list$X2 < 0.8] <- "intermediate"
  pred_res_level3_list$final[pred_res_level3_list$X2 < 0.5] <- "unassigned"
  
  rownames(pred_res_level3_list) <- unlist(lapply(strsplit(rownames(pred_res_level3_list), "\\."), function(x) {
    x <- paste(x[-1], collapse  = ".")
    x
  }))
  
  pred_res_level3 <- pred_res_level2
  
  pred_res_level3[rownames(pred_res_level3_list)] <- pred_res_level3_list$final
  
  print(table(pred_res_level3))
  
  pred_res = cbind(level1 = pred_res_level1,
                   level2 = pred_res_level2,
                   level3 = pred_res_level3)
  return(pred_res)
  
}




contstructTrainObjCV <- function(sce, celltype, upper = TRUE) {
  
  nms <- rownames(sce)
  if (upper) {
    bad_genes <- unique(c(grep("^MT-", nms, v=T), 
                          grep("^MTMR", nms, v=T), 
                          grep("^MTND", nms, v=T),
                          grep("RPL|RPS", nms, v=T),
                          "NEAT1","TMSB4X", "TMSB10"))
  } else {
    bad_genes <- unique(c(grep("^mt-", nms, v=T), 
                          grep("^Mtmr", nms, v=T), 
                          grep("^Mtnd", nms, v=T),
                          grep("Rpl|Rps", nms, v=T),
                          "Neat1","Tmsb4x", "Tmsb10"))
  }
  
  gene_corMat <- qlcMatrix::cosSparse(t((logcounts(sce)[bad_genes, ])), 
                                      t((logcounts(sce)[!rownames(sce) %in% bad_genes, ])))
  gene_corMat_max <- apply(gene_corMat, 2, max, na.rm = TRUE)
  exclude_genes <- c(bad_genes, names(gene_corMat_max)[gene_corMat_max > 0.7])
  
  
  
  
  set.seed(2020)
  library(cvTools)
  train_cv_list <- list()
  for (i in 1:10) {
    cat("=================== ")
    cat(i)
    cat("===================", "\n")
    
    cv_idx <- cvFolds(ncol(sce), K = 5)
    train_cv <- list()
    
    train_cv <- pbmclapply(1:5, function(k) {
      test_idx <- cv_idx$subsets[cv_idx$which == k]
      train_idx <- cv_idx$subsets[cv_idx$which != k]
      
      
      scClassify(exprsMat_train = as.matrix(logcounts(sce[, train_idx])[!rownames(sce) %in% exclude_genes, ]),
                 cellTypes_train = celltype[train_idx],
                 exprsMat_test = list(
                   test = as.matrix(logcounts(sce[, test_idx]))
                 ),
                 cellTypes_test = list(
                   test = celltype[test_idx]
                 ),
                 similarity = c("pearson"),
                 selectFeatures = c("limma"),
                 algorithm = c("WKNN"),
                 parallel = F,
                 prob_threshold = 0.7,
                 hopach_kmax = 5,
                 verbose = TRUE,
                 k = 10,
                 topN = 30)
      
      
    }, mc.cores = 5)
    
    train_cv_list[[i]] <- train_cv
  }
  
  pred_cv <- lapply(train_cv_list, function(y) unlist(lapply(y, function(x) 
    x$testRes$test$pearson_WKNN_limma$predRes))[colnames(sce)])
  
  pred_cv <- lapply(pred_cv, function(x) {
    
    pred_length <- unlist(lapply(strsplit(x, "_"), length))
    x[pred_length >= 2] <- "intermediate"
    x
    
  })
  
  pred_cv <- do.call(cbind, pred_cv)
  
  pred_cv_labels <- apply(pred_cv, 1, function(x) {
    tab <- table(x)/length(x)
    if (max(tab) >= 0.9) {
      names(tab)[which.max(tab)] 
    } else{
      "unassigned"
    }
    
    
  })
  
  print(mean(pred_cv_labels == celltype))
  
  train_idx <- pred_cv_labels == celltype
  
  sce_train <- sce[, train_idx]
  sce_train
  return(sce_train)
}


li <- readRDS("../data/melanoma_li.rds")
sadeFeldman <- readRDS("../data/melanoma_sadeFeldman.rds")

li
sadeFeldman

source("scripts/multiLevelscClassify.R")

dir.create("../results/")
dir.create("../results/melanoma/")
# Construct reference using Li data

if (file.exists("../results/melanoma/trainLi_scClassify.rds")) {
  trainLi_scClassify <- readRDS("../results/melanoma/trainLi_scClassify.rds")
} else {
  train_li_obj <- contstructTrainObjCV(li, li$cellTypes)
  
  tab <-  table(train_li_obj$cellTypes)
  keep_celltypes <- names(tab[tab >= 10])
  train_li_obj <- train_li_obj[, train_li_obj$cellTypes %in% keep_celltypes]
  trainLi_scClassify <- train_scClassify(exprsMat_train = logcounts(train_li_obj),
                                         cellTypes_train = train_li_obj$cellTypes)
  saveRDS(train_li_obj, file = "../results/melanoma/train_li_obj.rds")
  saveRDS(trainLi_scClassify, file = "../results/melanoma/trainLi_scClassify.rds")
}

if (file.exists("../results/melanoma/sadeFeldman_scClassify_res.rds")) {
  sadeFeldman_scClassify_res <- readRDS("../results/melanoma/sadeFeldman_scClassify_res.rds")
} else {
  sadeFeldman_scClassify_res <- runMultiLevelscClassify(sadeFeldman, trainLi_scClassify, 
                                                        n_cells_subset = 30000,
                                                        ncores = 1)
  saveRDS(sadeFeldman_scClassify_res, file = "../results/melanoma/sadeFeldman_scClassify_res.rds")
}


common_genes <- intersect(rownames(li), rownames(sadeFeldman))
exprsMat <- cbind(logcounts(li)[common_genes, ], logcounts(sadeFeldman)[common_genes, ])
celltype_label <- c(li$cellTypes, sadeFeldman$scClassify)
batch_label <- c(rep("li", ncol(li)),
                 rep("sadeFeldman", ncol(sadeFeldman)))

system.time(scMerge_res <- scMerge2(exprsMat, 
                                    batch_label, 
                                    celltype_label, SerialParam(),
                                    cosineNorm = FALSE, k_psuedoBulk = 5,
                                    pseudoBulk_fn = create_pseudoBulk_graph,
                                    ctl = rownames(exprsMat),
                                    ruvK = 30, ncores = 1))



hvg_fit <- scran::modelGeneVar(exprsMat, block = batch_label)
hvg <- scran::getTopHVGs(hvg_fit, n = 2000)

scMerge_res$newY[scMerge_res$newY < 0.0001] <- 0
scMerge_res$newY <- as(scMerge_res$newY, "dgCMatrix")
gc(reset = TRUE)


saveRDS(t(scMerge_res$newY[batch_label == "li", ]), 
        file = "../results/melanoma/exprsMat_scMerge2_li.rds")
saveRDS(t(scMerge_res$newY[batch_label != "li", ]), 
        file = "../results/melanoma/exprsMat_scMerge2_sadefeldman.rds")

meta_sadeFeldman <- colData(sadeFeldman)
saveRDS(meta_sadeFeldman, 
        file = "../results/melanoma/meta_sadefeldman.rds")

meta_li <- colData(li)
saveRDS(meta_li, 
        file = "../results/melanoma/meta_li.rds")

# For a new dataset, please modify the paths of gene expression matrix as well as the 
# celltype_info and patient_info
library(SingleCellExperiment)
library(CellChat)
library(dplyr)
future::plan("multiprocess", workers = 1)

# normalised gene expression matrix (gene by cell)
exprsMat <- readRDS("../results/melanoma/exprsMat_scMerge2_sadefeldman.rds")
# meta info for each cell
meta <- readRDS("../results/melanoma/meta_sadefeldman.rds")
# The cell type info used as cell grouping
celltype_info <- meta$scClassify
# The sample (patient) info used as cell grouping
patient_info <- meta$patient



sce <- SingleCellExperiment(assay = list(logcounts = exprsMat),
                            colData = meta)

sce$celltype <- celltype_info
sce$patient <- patient_info
tab <- table(sce$patient)
sample_name <- names(tab)[tab >= 100]

dir.create("../results/melanoma/CCI_sadefeldman/")

for (s in 1:length(sample_name)) {
  
  print(sample_name[s])
  
  sce_patient <- sce[, sce$patient == sample_name[s]]
  sce_patient <- sce_patient[, !sce_patient$celltype %in% 
                               c("intermediate", "unassigned")]
  
  
  tab <- table(sce_patient$celltype)
  print(tab)
  sce_patient <- sce_patient[, sce_patient$celltype %in% names(which(tab > 5))]
  print(table(sce_patient$celltype))
  
  if (length(table(sce_patient$celltype)) > 3) {

    
    cellchat <- createCellChat(object = as.matrix(logcounts(sce_patient)), 
                               meta = data.frame(colData(sce_patient)), 
                               group.by = "celltype")
    
    cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity

    
    levels(cellchat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
    
    CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
    showDatabaseCategory(CellChatDB)
    # Show the structure of the database
    dplyr::glimpse(CellChatDB$interaction)
    
    CellChatDB.use <- CellChatDB # use Secreted Signaling for cell-cell communication analysis
    cellchat@DB <- CellChatDB.use # set the used database in the object
    
    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
    # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    saveRDS(cellchat, file = paste("../results/melanoma/CCI_sadefeldman/cellchat_sadefeldman_", sample_name[s], ".rds", sep = ""))
    
  }
  
}










library(Seurat)
library(ggplot2)
library(NMF)
library(ggalluvial)
library(ggsci)
library(ComplexHeatmap)
library(CellChat)
library(future)
library(patchwork)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 4000 * 1024^2)
future::plan("multisession", workers = 4)

# Set working directory
setwd("/data/Blizard-AlazawiLab/rk/cellchat/BAFF")

# Prepare required input data for CellChat analysis
SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/BAFF/liver_seurat_HIGH_receptor.rds")

# Arrange the clusters in increasing order
Idents(SeuObj) <- "cluster"
SeuObj@active.ident <- factor (SeuObj@active.ident,
                                 levels = c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22','C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29'))


# Keep clusters having min 10 cells
SeuObjcellchat <- subset(x = SeuObj, idents = c("C23", "C24", "C26", "C29"), invert = TRUE)

data.input <- SeuObjcellchat[["RNA"]]$data # normalized data matrix
labels <- Idents(SeuObjcellchat)
# create a dataframe of the cell labels
meta <- data.frame(labels = labels, row.names = names(labels))

#### Create a CellChat object for liver_seurat_HIGH_receptor ####
cellchat <- createCellChat(object = SeuObjcellchat, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")

saveRDS(cellchat, file = "cellchat_liver_seurat_HIGH_receptor.rds")


#### Create a CellChat object for liver_seurat_LOW_receptor ####
# Prepare required input data for CellChat analysis
SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/BAFF/liver_seurat_LOW_receptor.rds")

# Arrange the clusters in increasing order
Idents(SeuObj) <- "cluster"
SeuObj@active.ident <- factor (SeuObj@active.ident,
                               levels = c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22','C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29'))


# Keep clusters having min 10 cells
SeuObjcellchat <- subset(x = SeuObj, idents = c("C23", "C24", "C26", "C29"), invert = TRUE)

data.input <- SeuObjcellchat[["RNA"]]$data # normalized data matrix
labels <- Idents(SeuObjcellchat)
# create a dataframe of the cell labels
meta <- data.frame(labels = labels, row.names = names(labels))

# Create a CellChat object for liver_seurat_LOW_receptor
cellchat <- createCellChat(object = SeuObjcellchat, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# Subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")

saveRDS(cellchat, file = "cellchat_liver_seurat_LOW_receptor.rds")
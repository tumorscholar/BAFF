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
future::plan("multisession", workers = 1)

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

# Compute probability per pathway
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate cell–cell communication
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "cellchat_liver_seurat_HIGH_receptor.rds")


#### Create a CellChat object for liver_seurat_LOW_receptor ####
# Prepare required input data for CellChat analysis
SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/BAFF/liver_seurat_LOW_receptor.rds")

# Arrange the clusters in increasing order
Idents(SeuObj) <- "cluster"
SeuObj@active.ident <- factor (SeuObj@active.ident,
                               levels = c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22','C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29'))


# Keep clusters having min 10 cells and common between groups
SeuObjcellchat <- subset(x = SeuObj, idents = c("C23", "C24", "C26", "C28", "C29"), invert = TRUE)

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

# Compute probability per pathway
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate cell–cell communication
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "cellchat_liver_seurat_LOW_receptor.rds")


#### load cellchat objects ####
cellchatHR <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/BAFF/cellchat_liver_seurat_HIGH_receptor.rds")
cellchatLR <- readRDS("/data/Blizard-AlazawiLab/rk/cellchat/BAFF/cellchat_liver_seurat_LOW_receptor.rds")

# Merge datasets for comparison
object.list <- list(BAFFreceptorHigh = cellchatHR, BAFFreceptorLow = cellchatLR)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Identify cell populations with significant changes in sending or receiving signals
num.link <- unlist(sapply(object.list, function(x) {
 rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
}))

weight.MinMax <- c(min(num.link, na.rm = TRUE), max(num.link, na.rm = TRUE)) 

gg <- list()
for (i in 1:length(object.list)) {
 gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                              title = names(object.list)[i], 
                                              weight.MinMax = weight.MinMax, 
                                              label.size = 6, font.size = 14, font.size.title = 14 ) +  
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) + 
  scale_x_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
  
  theme_classic(base_size = 12) +
  theme(
   axis.text.x = element_text(angle = 45, hjust = 1),
   axis.line = element_line(color = "black"),
   legend.position = "right",
   legend.direction = "vertical",
   legend.title = element_text(size = 10),
   legend.text = element_text(size = 8),
   legend.key.size = unit(1.5, "lines")
  )
}
ggp <- patchwork::wrap_plots(plots = gg)
ggp

ggsave("/data/home/hdx044/plots/BAFF/Differntial_interaction_strength_2D_liver.png", dpi = 300, width = 12, height = 6)

#### Compare the total number of interactions and interaction strength ####
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

ggp <- gg1 +
 theme_classic(base_size = 12) +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.line = element_line(color = "black"),
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  legend.key.size = unit(1.5, "lines")
 )

ggp

# Save high-resolution SVG
ggsave(
 filename = "cellcahtBAFFinteractions.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/BAFF",
 width = 3,     
 height = 4,
 units = "in"
)

# Identify dysfunctional signaling by using differential expression analysis

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "BAFFreceptorHigh"

# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)

# extract the ligand-receptor pairs with upregulated ligands in BAFFreceptorHigh
net.up <- subsetCommunication(cellchat, net = net, datasets = "BAFFreceptorHigh",ligand.logFC = 0.05, receptor.logFC = NULL)

write.csv(net.up, '/data/home/hdx044/files/cellchat/BAFF/net_up_baff_receptor_high.csv')

# extract the ligand-receptor pairs with downregulated ligands in BAFFreceptorHigh
net.down <- subsetCommunication(cellchat, net = net, datasets = "BAFFreceptorHigh",ligand.logFC = -0.05, receptor.logFC = NULL)

write.csv(net.down, '/data/home/hdx044/files/cellchat/BAFF/net_down_baff_receptor_high.csv')

# Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs
pairLR.use.up = net.up[, "interaction_name", drop = F]

# Comparing communications on a merged object
gg1 <- netVisual_bubble(cellchat, sources.use = c(9), targets.use = c(1:24),  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signalling", angle.x = 90, remove.isolate = T, dot.size.min = 4,
                        dot.size.max = 5,
                        font.size = 16,
                        font.size.title = 14,
                        color.text = c("maroon", "mistyrose"))
# Increase legend font size
gg1 <- gg1 + theme(
 legend.title = element_text(size = 16),   # title font
 legend.text  = element_text(size = 14),   # text font
 legend.key.size = unit(1.5, "lines")      # enlarge legend symbols
)

gg1

# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)








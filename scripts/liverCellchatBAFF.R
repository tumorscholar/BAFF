library(Seurat)
library(ggplot2)
library(NMF)
library(ggalluvial)
library(ggsci)
library(ComplexHeatmap)
library(CellChat)
library(future)
library(patchwork)
library(officer)
library(rvg)
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

#### Identify dysfunctional signaling by using differential expression analysis ####

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

#### Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs ####

# Check pval distribution
cat("pval summary:\n")
print(summary(net.up$pval))

cat("\nHow many interactions with pval = 0:\n")
print(sum(net.up$pval == 0))

cat("\nTotal interactions:\n")
print(nrow(net.up))

cat("\npval = 0 percentage:\n")
print(sum(net.up$pval == 0) / nrow(net.up) * 100)

# Check if any non-zero pvals exist
cat("\nNon-zero pvals:\n")
print(net.up[net.up$pval > 0, c("source", "target", "interaction_name_2", "pathway_name", "prob", "pval")])

# Check fibroblast-specific interactions without prob threshold
cat("All interactions targeting Fibroblasts (C13):\n")
print(net.up[as.character(net.up$target) == "C13" & 
              as.character(net.up$source) %in% sources_of_interest,
             c("source", "target", "interaction_name_2", "pathway_name", "prob", "pval")])



# Cluster to cell type mapping
cluster_to_type <- c(
 "C0"  = "CD8 Effector Memory",
 "C1"  = "Naïve-Activated T cells",
 "C2"  = "CD56 dim NK cells",
 "C3"  = "CD14 Monocytes",
 "C4"  = "CD8 Cytotoxic T cells",
 "C5"  = "Macrophage",
 "C6"  = "Blood vessels",
 "C7"  = "CD4 T cells",
 "C8"  = "Naïve B cells",
 "C9"  = "ACKR1 Endothelial cells",
 "C10" = "Myeloid Dendritic cells",
 "C11" = "CD16 Monocytes",
 "C12" = "CD56 high NK cells",
 "C13" = "Fibroblasts",
 "C14" = "Plasmacytoid Dendritic cells",
 "C15" = "Mast cells",
 "C16" = "TAGLN Endothelial cells",
 "C17" = "Erythroid cells",
 "C18" = "Liver Ductal cells",
 "C19" = "Plasma cells_1",
 "C20" = "Hepatocytes",
 "C21" = "CD4 Proliferating T cells",
 "C22" = "CD4 Dendritic cells",
 "C23" = "TREM2 Dendritic cells",
 "C24" = "Mesothelial cells",
 "C25" = "Adipocytes_1",
 "C26" = "Adipocytes_2",
 "C27" = "B cells",
 "C28" = "Platelets",
 "C29" = "Plasma cells_2"
)

# Source labels
source_labels <- c(
 "C8"  = "Naïve B cells (C8)",
 "C19" = "Plasma cells_1 (C19)",
 "C22" = "CD4 Dendritic cells (C22)",
 "C27" = "B cells (C27)"
)

pathways_of_interest <- c(
 # BAFF family - core to your story
 "BAFF",          # BAFF ligand - B cell survival
 "APRIL",         # APRIL - BAFF family, plasma cell survival
 
 # Pro-fibrotic growth factors
 "BMP",           # BMP6 - HSC activation, iron metabolism
 "WNT",           # WNT10A - fibrosis/HSC
 "ncWNT",         # WNT5B - non-canonical fibrosis
 "IGF",           # IGF1 - hepatocyte survival/fibrosis
 "VEGF",          # VEGFB - angiogenesis/fibrosis
 "MK",            # Midkine - HSC activation
 "GRN",           # Granulin - liver fibrosis driver
 "GAS",           # GAS6 - TAM receptor, fibrosis resolution
 
 # Inflammatory mediators
 "MIF",           # MIF - macrophage migration, inflammation
 "IL10",          # IL10 - immunoregulation
 "IL16",          # IL16 - T cell chemoattractant
 "CCL",           # CCL3 - immune recruitment
 "CXCL",          # CXCL16 - NK/T cell recruitment
 "TWEAK",         # TNF family - hepatocyte death/fibrosis
 
 # Cell adhesion / matrix
 "ICAM",          # ICAM2 - leukocyte adhesion
 "VCAM",          # VCAM - immune cell recruitment
 "PECAM1",        # PECAM1 - endothelial/immune adhesion
 "PECAM2",        # PECAM1-CD38 contact
 "NECTIN",        # Nectin - cell adhesion
 "SELL",          # Selectin - leukocyte rolling/recruitment
 "ADGRE",         # ADGRE5 - immune cell contact
 "ADGRA",         # ADGRA2 - SDC1 receptor, plasma cell niche
 
 # Immune regulation
 "MHC-I",         # Antigen presentation to CD8 T cells
 "MHC-II",        # Antigen presentation to CD4 T cells
 "LAIR1",         # LAIR1 - inhibitory receptor
 "SEMA4",         # Semaphorin - immune regulation
 "SEMA7",         # Semaphorin - fibrosis associated
 
 # Signalling
 "CypA",          # Cyclophilin A - inflammation
 "Prostaglandin", # PGE2 - inflammation/immune modulation
 "CD99",          # CD99 - leukocyte migration
 "GAS"            # GAS6 - TAM receptor (duplicate removed below)
)

# Remove duplicates
pathways_of_interest <- unique(pathways_of_interest)

# Remove prob threshold - filter by pathway and source only
net_filtered <- net.up %>%
 filter(as.character(source) %in% sources_of_interest,
        pathway_name %in% pathways_of_interest) %>%
 mutate(
  source_label = source_labels[as.character(source)],
  target_label = cluster_to_type[as.character(target)]
 ) %>%
 group_by(source, interaction_name) %>%
 mutate(max_prob = max(prob)) %>%
 ungroup()

# Sanity check
cat("Fibroblast interactions:\n")
print(net_filtered[as.character(net_filtered$target) == "C13",
                   c("source", "interaction_name_2", "pathway_name", "prob")])

cat("\nTotal interactions per source:\n")
print(table(net_filtered$source_label, useNA = "ifany"))

# Build bubble data
bubble_data_target <- net_filtered %>%
 filter(!is.na(target_label)) %>%
 group_by(source_label, target_label, interaction_name_2, pathway_name) %>%
 summarise(
  prob       = mean(prob),
  ligand_pct = mean(ligand.pct.1, na.rm = TRUE),
  .groups    = "drop"
 ) %>%
 mutate(source_label = factor(source_label,
                              levels = c("Naïve B cells (C8)",
                                         "Plasma cells_1 (C19)",
                                         "CD4 Dendritic cells (C22)",
                                         "B cells (C27)")))

# Plot
# Split by source and plot individually
source_levels <- c("Naïve B cells (C8)", "Plasma cells_1 (C19)",
                   "CD4 Dendritic cells (C22)", "B cells (C27)")

# Get global prob range across all sources for consistent scale
prob_range <- range(bubble_data_target$prob, na.rm = TRUE)
pct_range  <- range(bubble_data_target$ligand_pct, na.rm = TRUE)

plot_list <- list()

for (src in source_levels) {
 
 dat <- bubble_data_target %>% filter(source_label == src)
 
 p <- ggplot(dat,
             aes(x = interaction_name_2,
                 y = target_label,
                 size = ligand_pct,
                 color = prob)) +
  geom_vline(aes(xintercept = interaction_name_2),
             color = "grey85", linewidth = 0.3) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(range  = c(1, 6),
                        limits = pct_range,    # consistent size scale
                        name   = "Ligand expression\n(% source cells)") +
  scale_color_viridis_c(
   option    = "magma",
   direction = 1,                             # force consistent direction
   limits    = prob_range,                    # consistent colour scale
   name      = "Communication\nprobability"
  ) +
  facet_grid(. ~ pathway_name, scales = "free_x", space = "free_x") +
  labs(
   title = paste("Upregulated interactions —", src),
   x     = "Ligand - Receptor",
   y     = "Target cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
   axis.text.x        = element_text(angle = 90, hjust = 1, 
                                     vjust = 0.5, size = 14),
   axis.text.y        = element_text(size = 10),
   axis.title         = element_text(size = 14, face = "bold"),
   strip.background   = element_rect(fill = "grey90", color = NA),
   strip.text.x       = element_text(face = "bold", size = 9, angle = 90,
                                     hjust = 0, vjust = 0.5),
   plot.title         = element_text(hjust = 0.5, face = "bold", size = 14),
   panel.spacing      = unit(0.3, "lines"),
   legend.position    = "right",
   legend.title       = element_text(size = 14),
   legend.text        = element_text(size = 14),
   panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
   panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )
 
 plot_list[[src]] <- p
 
 # Save individual plots
 ggsave(
  filename = file.path(plot_dir, paste0("cellchat_bubble_",
                                        gsub(" ", "_", gsub("[^A-Za-z0-9 ]", "", src)),
                                        ".png")),
  plot   = p,
  width  = 16,
  height = 8,
  units  = "in",
  dpi    = 300
 )
 
 cat("Saved:", src, "\n")
}

#### Pathway plot ####


# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

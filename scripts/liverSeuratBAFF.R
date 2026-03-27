# Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(SeuratExtend)
library(escape)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(viridisLite)
library(ggrepel)

# Load seurat object
SeuObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

DefaultAssay(SeuObj) <- "RNA"

baff_signaling_genes <- c(
 "TNFRSF13C",  # BAFF-R (BAFF receptor)
 "TNFRSF13B",  # TACI
 "TNFRSF17",   # BCMA
 "TNFSF13B"    # BAFF (B cell–activating factor)
)

# Check Baff expression in Tissue
DotPlot2(SeuObj, features = "TNFSF13B", split.by = 'Tissue')


# Arrange Stage in the desired order
SeuObj$Stage <- factor(
 SeuObj$Stage,
 levels = c("Healthy", "No_fibrosis", "Fibrosis")
)

ggp <- DotPlot2(SeuObj, features = baff_signaling_genes, color_scheme = "BuRd", group.by = "Tissue", split.by = "Stage") +
 ggtitle("BAFF and its receptors") +
 coord_flip()

ggp 

# Save high-resolution SVG
ggsave(
 filename = "BAFFandReceptorsTissue.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/BAFF",
 width = 5,     
 height = 5,
 units = "in"
)

# Subset LIVER
LIVER <- subset(SeuObj, subset= Tissue %in% "LIVER")

baff_signaling_genes <- c(
 "TNFRSF13C",  # BAFF-R (BAFF receptor)
 "TNFRSF13B",  # TACI (Transmembrane activator and CAML interactor)
 "TNFRSF17",   # BCMA
 "TNFSF13B"    # BAFF (B cell–activating factor)
 )

LIVER$cell_type_with_cluster

# Extract all labels
labs <- LIVER$cell_type_with_cluster

# Extract cluster numbers
cluster_num <- as.numeric(sub(".*\\(C([0-9]+)\\).*", "\\1", labs))

# Build a data frame of unique labels + cluster numbers
df_unique <- unique(data.frame(
 label = labs,
 cluster = cluster_num
))

# Sort unique labels by cluster number
df_unique <- df_unique[order(df_unique$cluster), ]

# Create sorted factor
labs_sorted <- factor(labs, levels = df_unique$label)

# Assign back to Seurat metadata
LIVER$cell_type_with_cluster <- labs_sorted
Idents(LIVER) <- "cell_type_with_cluster"
levels(LIVER)

# Ploting UMAP for cell_type_with_cluster
ggp <- DimPlot(LIVER, reduction = "umap", raster = FALSE, label = F, group.by = "cell_type_with_cluster") +
 ggtitle("Liver clusters") +
 theme_classic(base_size = 12) +
 theme(
  axis.text.x = element_text(angle = 90, hjust = 1),
  axis.line = element_line(color = "black"),
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.key.size = unit(0.8, "lines"),  # reduce point size in legend
  legend.spacing.y = unit(0.8, "cm")     # reduce spacing between items
 ) +
 guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))  # one column legend, small dots

ggp

# Final TIFF for journal
setwd("/data/home/hdx044/plots/BAFF")
tiff("liver_cell_type.tiff", width=22, height=15, units="cm", res=600, compression="lzw")
print(ggp)
dev.off()


ggp <- DotPlot2(LIVER, features = baff_signaling_genes, color_scheme = "BuRd") +
 ggtitle("BAFF and its receptors") +
 coord_flip()+ 
 guides(
  color = guide_colorbar(order = 1),  # Average Expression on top
  size = guide_legend(order = 2))      # Percent Expressed below

ggp 

# Save high-resolution SVG
ggsave(
 filename = "BAFFandReceptors.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/BAFF",
 width = 5,     
 height = 5,
 units = "in"
)


baff_signaling_ADTs <- c(
 "Hu.CD268",  #  TNFRSF13C also BAFF-R (BAFF receptor)
 "Hu.CD267"  #  TNFRSF17 also BCMA (B cell maturation antigen)
)

DefaultAssay(LIVER) <- "ADTonly"

ggp <- DotPlot2(LIVER, features = baff_signaling_ADTs, color_scheme = "BuRd") +
 ggtitle("BAFF and its receptors") +
 coord_flip()+ 
 guides(
  color = guide_colorbar(order = 1),  # Average Expression on top
  size = guide_legend(order = 2))      # Percent Expressed below

ggp 

# Save high-resolution SVG
ggsave(
 filename = "BAFFReceptorsADTs.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/BAFF",
 width = 5,     
 height = 5,
 units = "in"
)

#### Generate count plots Stage wise ####
# Set identities
Idents(LIVER) <- 'cluster'

#  Order clusters in metadata column directly
LIVER$cluster <- factor(
 LIVER$cluster,
 levels = c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10',
            'C11','C12','C13','C14','C15','C16','C17','C18','C19','C20',
            'C21','C22','C23','C24','C25','C26','C27','C28','C29')
)

# ALSO order the active.ident (for downstream Seurat functions)
LIVER@active.ident <- factor(
 LIVER@active.ident,
 levels = c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10',
            'C11','C12','C13','C14','C15','C16','C17','C18','C19','C20',
            'C21','C22','C23','C24','C25','C26','C27','C28','C29')
)

# Order Stage
LIVER$Stage <- factor(
 LIVER$Stage,
 levels = c("Healthy", "No_fibrosis", "Fibrosis")
)

# Create summary - clusters will now be ordered
summary <- LIVER@meta.data %>%
 group_by(cluster, Stage) %>%
 summarise(count = n(), .groups = 'drop') %>%
 group_by(cluster) %>%
 mutate(proportion = count / sum(count)) %>%
 ungroup()

# Create color palette
viridis_cols <- viridis(length(unique(summary$Stage)))
names(viridis_cols) <- levels(LIVER$Stage)  # Use levels, not unique

# Plot
ggp <- ggplot(summary, aes(x = cluster, y = proportion, fill = Stage)) +
 geom_bar(stat = "identity", position = "stack") +
 scale_fill_manual(values = viridis_cols) +
 scale_y_continuous(expand = c(0, 0)) +
 labs(
  title = "Proportions of cells in each cluster by stage",
  x = "Cluster", 
  y = "Proportion of Cells"
 ) +
 theme_classic(base_size = 12) +
 theme(
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  axis.line = element_line(color = "black"),
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.key.size = unit(0.8, "lines"),
  legend.spacing.y = unit(0.5, "cm")
 ) +
 guides(fill = guide_legend(override.aes = list(size = 3), ncol = 1))

# Display plot
print(ggp)

setwd("/data/home/hdx044/plots/BAFF")
tiff("liver_cluster_Proportion_by_Stage.tiff", width=25, height=10, units="cm", res=600, compression="lzw")
print(ggp)
dev.off()


# Save high-resolution SVG
ggsave(
 filename = "liver_cluster_Proportion_by_Stage.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/BAFF",
 width = 5,     
 height = 3,
 units = "in",
 dpi = 600
)


saveRDS(LIVER, '/data/Blizard-AlazawiLab/rk/seurat/LIVER.rds')

# BAFF and its receptors clusters
liverBaff <- subset(LIVER, subset = cluster %in% c("C3", "C5", "C8", "C10", "C11", "C13", "C14", "C19", "C22", "C23", "C27", "C29"))

saveRDS(liverBaff, '/data/Blizard-AlazawiLab/rk/seurat/liverBaff.rds')

#### Load liver BAFF obj ####
liverBaff <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/liverBaff.rds')

Idents(liverBaff) <- "StageSep"
DefaultAssay(liverBaff) <- "RNA"

# Arrange Stage in the desired order
liverBaff$StageSep <- factor(
 liverBaff$StageSep,
 levels = c("H", "F0", "F1", "F2", "F3")
)

ggp <- DotPlot2(liverBaff, features = baff_signaling_genes, color_scheme = "BuRd", group.by = "StageSep") +
 ggtitle("BAFF and its receptors") +
 guides(
  color = guide_colorbar(order = 1),  # Average Expression on top
  size = guide_legend(order = 2))      # Percent Expressed below

ggp 

# Save high-resolution SVG
ggsave(
 filename = "LiverBAFFandReceptorsStageSep.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/BAFF",
 width = 5,     
 height = 5,
 units = "in"
)


baff_signaling_ADTs <- c(
 "Hu.CD268",  #  TNFRSF13C also BAFF-R (BAFF receptor)
 "Hu.CD267"  #  TNFRSF17 also BCMA (B cell maturation antigen)
)

DefaultAssay(liverBaff) <- "ADTonly"

ggp <- DotPlot2(liverBaff, features = baff_signaling_ADTs, color_scheme = "BuRd", group.by = "StageSep") +
 ggtitle("BAFF and its receptors") +
 guides(
  color = guide_colorbar(order = 1),  # Average Expression on top
  size = guide_legend(order = 2))      # Percent Expressed below

ggp 

# Save high-resolution SVG
ggsave(
 filename = "LiverBAFFandReceptorsADTstage.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/BAFF",
 width = 3,     
 height = 4,
 units = "in"
)

#### Make BAFF expressing and absent cell obj ####

DefaultAssay(liverBaff) <- "RNA"
# Keep only genes present
genes_present <- intersect("TNFSF13B", rownames(liverBaff))
if (length(genes_present) == 0) {
 stop("BAFF gene is not present in the object.")
}

# Get expression matrix (assumes Seurat and using normalized 'data' slot; adjust if needed)
expr <- tryCatch(
 GetAssayData(liverBaff, slot = "data")[genes_present, , drop = FALSE],
 error = function(e) liverBaff[genes_present, , drop = FALSE]  # if it's a plain matrix/dgCMatrix
)

# Compute per-cell max expression across the BAFF receptor genes
baff_expr <- matrixStats::colMaxs(as.matrix(expr))

# Label High if any receptor >= 0.25, otherwise Low
liverBaff$BAFFstatus <- ifelse(baff_expr >= 0.25, "High", "Low")

# Quick sanity checks
table(liverBaff$BAFFstatus)
# High  Low 
# 3371 7745


# check obj
Idents(liverBaff) <- 'BAFFstatus'
VlnPlot2(liverBaff, features = "TNFSF13B")

#### Volcano plot HIGH vs LOW BAFF expressing cells Cells ####

# Step 1: Differential Expression Analysis
# Set identity to BAFFreceptorStatus
Idents(liverBaff) <- "BAFFstatus"

# Perform DE analysis (HIGH vs LOW)
DE_results <- FindMarkers(
 liverBaff,
 ident.1 = "High",
 ident.2 = "Low",
 test.use = "wilcox",  # Wilcoxon rank-sum test
 logfc.threshold = 0,  # Include all genes for volcano
 min.pct = 0,          # Include all genes
 verbose = TRUE
)

# Add gene names as column
DE_results$gene <- rownames(DE_results)

# Step 2: Add Classification Column

DE_results <- DE_results %>%
 mutate(
  direction = case_when(
   avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Upregulated in HIGH",
   avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Upregulated in LOW",
   TRUE ~ "Not Significant"
  ),
  neg_log10_padj = -log10(p_val_adj)
 )

# Export significant genes only
DE_sig <- DE_results %>%
 filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
 arrange(desc(abs(avg_log2FC)))

write.csv(DE_sig, "/data/home/hdx044/files/BAFF/DE_HIGH_vs_LOW_BAFF_significant.csv")


# Genes to label on volcano plot

label_gene_names <- c(
 # BAFF receptors (confirmation of stratification)
 "TNFRSF13B",  # TACI - BAFF receptor, confirms HIGH stratification
 "TNFRSF17",   # BCMA - BAFF receptor, confirms HIGH stratification
 "TNFRSF13C",  # BAFF-R - BAFF receptor, confirms HIGH stratification
 
 # Fibrosis-related - Upregulated in HIGH
 "BMP6",       # Iron metabolism + HSC activation in MASH
 "TGFBR3L",    # TGF-beta receptor signalling
 "WNT5B",      # Wnt signalling + fibrosis
 "IGF1",       # Hepatocyte survival + fibrosis
 "CTHRC1",     # ECM remodelling + HSC activation
 "SDC1",       # Plasma cell marker + fibrosis niche (CD138)
 "COL4A3",     # Collagen ECM component
 
 # Fibrosis-related - Upregulated in LOW
 "TGFB1",      # Canonical fibrosis driver
 "TIMP3",      # ECM regulation + TGF-beta pathway
 "CXCL12"      # Fibrosis niche + B cell homing
)

label_genes <- DE_results %>%
 filter(gene %in% label_gene_names)

print(table(DE_results$direction))

# Step 4: Create Volcano Plot
colors <- c(
 "Upregulated in HIGH" = "#E41A1C",
 "Upregulated in LOW"  = "#377EB8",
 "Not Significant"     = "grey70"
)

p_volcano <- ggplot(DE_results, aes(x = avg_log2FC, y = neg_log10_padj)) +
 
 geom_point(aes(color = direction), alpha = 0.6, size = 1.5) +
 
 geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
            color = "black", size = 0.5) +
 geom_hline(yintercept = -log10(0.05), linetype = "dashed",
            color = "black", size = 0.5) +
 
 geom_text_repel(
  data = label_genes,
  aes(label = gene, color = direction),
  size = 3,
  max.overlaps = 20,
  box.padding = 0.5,
  point.padding = 0.3,
  segment.color = "grey50",
  segment.size = 0.3,
  min.segment.length = 0,
  show.legend = FALSE,
  bg.color = "white",      # white box behind text
  bg.r = 0.15              # border radius of box
 ) +
 scale_color_manual(values = colors) +
 
 theme_classic(base_size = 12) +
 theme(
  legend.position = "top",
  legend.title = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
  axis.line = element_line(color = "black"),
  panel.grid.major = element_line(color = "grey90", size = 0.2)
 ) +
 
 labs(
  title = "Differential gene: High vs Low expressing BAFF receptors B cells",
  x = "log2 Fold Change (High vs Low)",
  y = "-log10(Adjusted P-value)",
  caption = "Dashed lines: |logFC| > 0.5, adjusted p < 0.05"
 ) +
 
 xlim(c(min(DE_results$avg_log2FC) - 0.5, max(DE_results$avg_log2FC) + 0.5)) +
 ylim(c(0, max(DE_results$neg_log10_padj) * 1.1))

print(p_volcano)

# Save plot
ggsave(
 "/data/home/hdx044/plots/BAFF/volcano_plot_HIGH_vs_LOW_BAFFreceptor.png", 
 p_volcano, 
 width = 12, 
 height = 5,
 dpi = 300
)

# Step 5: Export Results
# Export significant genes only
DE_sig <- DE_results %>%
 filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
 arrange(desc(abs(avg_log2FC)))

write.csv(DE_sig, "/data/home/hdx044/files/BAFF/DE_HIGH_vs_LOW_BAFFreceptor_significant.csv", row.names = FALSE)
cat("Saved: DE_HIGH_vs_LOW_BAFFreceptor_significant.csv\n")

#### Pathway analysis HIGH vs LOW BAFF Cells ####
# Getting genes set
C2 <- getGeneSets(library = "C2")

# Filter gene sets with names related to Reactome (case-sensitive)
reactome_sets <- C2[sapply(names(C2), function(x) grepl("REACTOME", x))]

# Check the filtered Reactome sets
length(reactome_sets)  # Number of Reactome gene sets

DefaultAssay(liverBaff) <- "RNA"

# Run Escape on Reactome pathway
liverBaff <- runEscape(liverBaff,
                                  method = "UCell",
                                  gene.sets = reactome_sets, 
                                  groups = 5000, 
                                  min.size = 0,
                                  new.assay.name = "escape.UCell")

saveRDS(liverBaff, '/data/Blizard-AlazawiLab/rk/seurat/liverBaff.rds')

ucell_mat <- GetAssayData(
 liverBaff,
 assay = "escape.UCell",
 slot = "data"
)

ucell_df <- as.data.frame(t(ucell_mat))
ucell_df$cell_id <- rownames(ucell_df)

# Add celltype to ucell_df
ucell_df$cell_type <- liverBaff$cell_type

# Calculate mean pathway scores by BAFF receptor status
pathway_BAFFmeans <- ucell_df %>%
 group_by(cell_type) %>%
 summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# Check the result
print(pathway_BAFFmeans)

# Pivot to long format - NOTE: No 'cluster' column anymore!
long_means <- pathway_BAFFmeans %>%
 pivot_longer(
  -cell_type,
  names_to = "pathway",
  values_to = "mean_score"
 )

# Check the result
print(head(long_means, 20))


#### Create BAFF receptor expressing and absent cell obj ####
# BAFF and its receptors clusters

liverBaffreceptorObj <- subset(liverBaff , subset = seurat_clusters %in% c('8', '19', '22', '27', '29'))

DefaultAssay(liverBaffreceptorObj) <- "RNA"

baff_receptors <- c("TNFRSF13C", "TNFRSF13B", "TNFRSF17")

# Keep only genes present
genes_present <- intersect(baff_receptors, rownames(liverBaffreceptorObj))
if (length(genes_present) == 0) {
 stop("None of the BAFF receptor genes are present in the object.")
}

# Get expression matrix (assumes Seurat and using normalized 'data' slot; adjust if needed)
expr <- tryCatch(
 GetAssayData(liverBaffreceptorObj, slot = "data")[genes_present, , drop = FALSE],
 error = function(e) liverBaffreceptorObj[genes_present, , drop = FALSE]  # if it's a plain matrix/dgCMatrix
)

# Compute per-cell max expression across the BAFF receptor genes
baff_expr <- matrixStats::colMaxs(as.matrix(expr))

# Label High if any receptor >= 0.25, otherwise Low
liverBaffreceptorObj$BAFFreceptorStatus <- ifelse(baff_expr >= 0.25, "High", "Low")

# Quick sanity checks
table(liverBaffreceptorObj$BAFFreceptorStatus)

# check obj
Idents(liverBaffreceptorObj) <- 'BAFFreceptorStatus'
VlnPlot2(liverBaffreceptorObj, features = baff_receptors)


#### Volcano plot HIGH vs LOW BAFF Receptor B Cells ####

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

cat("Performing differential expression analysis...\n")

# Step 1: Differential Expression Analysis
# Set identity to BAFFreceptorStatus
Idents(liverBaffreceptorObj) <- "BAFFreceptorStatus"

# Perform DE analysis (HIGH vs LOW)
DE_results <- FindMarkers(
 liverBaffreceptorObj,
 ident.1 = "High",
 ident.2 = "Low",
 test.use = "wilcox",  # Wilcoxon rank-sum test
 logfc.threshold = 0,  # Include all genes for volcano
 min.pct = 0,          # Include all genes
 verbose = TRUE
)

# Add gene names as column
DE_results$gene <- rownames(DE_results)

# Step 2: Add Classification Column

DE_results <- DE_results %>%
 mutate(
  direction = case_when(
   avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Upregulated in HIGH",
   avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Upregulated in LOW",
   TRUE ~ "Not Significant"
  ),
  neg_log10_padj = -log10(p_val_adj)
 )

# Summary
cat("\nDifferential Expression Summary:\n")

# Genes to label on volcano plot

label_gene_names <- c(
 # BAFF receptors (confirmation of stratification)
 "TNFRSF13B",  # TACI - BAFF receptor, confirms HIGH stratification
 "TNFRSF17",   # BCMA - BAFF receptor, confirms HIGH stratification
 "TNFRSF13C",  # BAFF-R - BAFF receptor, confirms HIGH stratification
 
 # Fibrosis-related - Upregulated in HIGH
 "BMP6",       # Iron metabolism + HSC activation in MASH
 "TGFBR3L",    # TGF-beta receptor signalling
 "WNT5B",      # Wnt signalling + fibrosis
 "IGF1",       # Hepatocyte survival + fibrosis
 "CTHRC1",     # ECM remodelling + HSC activation
 "SDC1",       # Plasma cell marker + fibrosis niche (CD138)
 "COL4A3",     # Collagen ECM component
 
 # Fibrosis-related - Upregulated in LOW
 "TGFB1",      # Canonical fibrosis driver
 "TIMP3",      # ECM regulation + TGF-beta pathway
 "CXCL12"      # Fibrosis niche + B cell homing
)

label_genes <- DE_results %>%
 filter(gene %in% label_gene_names)

print(table(DE_results$direction))

# Step 4: Create Volcano Plot
colors <- c(
 "Upregulated in HIGH" = "#E41A1C",
 "Upregulated in LOW"  = "#377EB8",
 "Not Significant"     = "grey70"
)

p_volcano <- ggplot(DE_results, aes(x = avg_log2FC, y = neg_log10_padj)) +
 
 geom_point(aes(color = direction), alpha = 0.6, size = 1.5) +
 
 geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
            color = "black", size = 0.5) +
 geom_hline(yintercept = -log10(0.05), linetype = "dashed",
            color = "black", size = 0.5) +
 
 geom_text_repel(
  data = label_genes,
  aes(label = gene, color = direction),
  size = 3,
  max.overlaps = 20,
  box.padding = 0.5,
  point.padding = 0.3,
  segment.color = "grey50",
  segment.size = 0.3,
  min.segment.length = 0,
  show.legend = FALSE,
  bg.color = "white",      # white box behind text
  bg.r = 0.15              # border radius of box
 ) +
 scale_color_manual(values = colors) +
 
 theme_classic(base_size = 12) +
 theme(
  legend.position = "top",
  legend.title = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
  axis.line = element_line(color = "black"),
  panel.grid.major = element_line(color = "grey90", size = 0.2)
 ) +
 
 labs(
  title = "Differential gene: High vs Low expressing BAFF receptors B cells",
  x = "log2 Fold Change (High vs Low)",
  y = "-log10(Adjusted P-value)",
  caption = "Dashed lines: |logFC| > 0.5, adjusted p < 0.05"
 ) +
 
 xlim(c(min(DE_results$avg_log2FC) - 0.5, max(DE_results$avg_log2FC) + 0.5)) +
 ylim(c(0, max(DE_results$neg_log10_padj) * 1.1))

print(p_volcano)

# Save plot
ggsave(
 "/data/home/hdx044/plots/BAFF/volcano_plot_HIGH_vs_LOW_BAFFreceptor.png", 
 p_volcano, 
 width = 12, 
 height = 5,
 dpi = 300
)

# Step 5: Export Results
# Export significant genes only
DE_sig <- DE_results %>%
 filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
 arrange(desc(abs(avg_log2FC)))

write.csv(DE_sig, "/data/home/hdx044/files/BAFF/DE_HIGH_vs_LOW_BAFFreceptor_significant.csv", row.names = FALSE)
cat("Saved: DE_HIGH_vs_LOW_BAFFreceptor_significant.csv\n")


#### Pathway analysis HIGH vs LOW BAFF Receptor B Cells ####
# Getting genes set
C2 <- getGeneSets(library = "C2")

# Filter gene sets with names related to Reactome (case-sensitive)
reactome_sets <- C2[sapply(names(C2), function(x) grepl("REACTOME", x))]

# Check the filtered Reactome sets
length(reactome_sets)  # Number of Reactome gene sets

# Run Escape on Reactome pathway
liverBaffreceptorObj <- runEscape(liverBaffreceptorObj,
                                  method = "UCell",
                                  gene.sets = reactome_sets, 
                                  groups = 5000, 
                                  min.size = 0,
                                  new.assay.name = "escape.UCell")

saveRDS(liverBaffreceptorObj, '/data/Blizard-AlazawiLab/rk/seurat/liverBaffreceptorObj.rds')

#### Load liver BAFF receptor obj ####
liverBaffreceptorObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/liverBaffreceptorObj.rds')

ucell_mat <- GetAssayData(
 liverBaffreceptorObj,
 assay = "escape.UCell",
 slot = "data"
)

ucell_df <- as.data.frame(t(ucell_mat))
ucell_df$cell_id <- rownames(ucell_df)

# Add BAFFreceptorStatus to ucell_df
ucell_df$BAFFreceptorStatus <- liverBaffreceptorObj$BAFFreceptorStatus

# Calculate mean pathway scores by BAFF receptor status
pathway_BAFFreceptorStatus_means <- ucell_df %>%
 group_by(BAFFreceptorStatus) %>%
 summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# Check the result
print(pathway_BAFFreceptorStatus_means)

# Pivot to long format - NOTE: No 'cluster' column anymore!
long_means <- pathway_BAFFreceptorStatus_means %>%
 pivot_longer(
  -BAFFreceptorStatus,  # CHANGED from -cluster
  names_to = "pathway",
  values_to = "mean_score"
 )

# Check the result
print(head(long_means, 20))

# Define UCell Cutoffs

UCELL_STRONG <- 0.20    # Strong enrichment
UCELL_MODERATE <- 0.15  # Moderate enrichment  
UCELL_WEAK <- 0.10      # Weak enrichment

cat("UCell Cutoffs (based on literature):\n")
cat("  Strong:   > 0.20\n")
cat("  Moderate: > 0.15\n")
cat("  Weak:     > 0.10\n\n")

# Get Max Scores Per Pathway

pathway_max_scores <- long_means %>%
 group_by(pathway) %>%
 summarise(
  max_score = max(mean_score),
  mean_score = mean(mean_score)
 ) %>%
 arrange(desc(max_score))

# Filter Meaningful Pathways

# Keep pathways that reach moderate threshold in at least one cluster
pathways_meaningful <- pathway_max_scores %>%
 filter(max_score > UCELL_MODERATE) %>%
 pull(pathway)

cat(sprintf("Pathways > 0.15: %d\n\n", length(pathways_meaningful)))

#### UCell pathways heatmap ####

ucell_mat <- GetAssayData(
 liverBaffreceptorObj,
 assay = "escape.UCell",
 slot = "data"
)

ucell_df <- as.data.frame(t(ucell_mat))
ucell_df$BAFFreceptorStatus <- liverBaffreceptorObj$BAFFreceptorStatus[rownames(ucell_df)]

# Pathways exactly as you specified
# Define BAFF/fibrosis-relevant pathways and categories ---

pathway_groups <- list(
 
 "NF-kB / Inflammatory Signaling" = c(
  "REACTOME-NF-KB-IS-ACTIVATED-AND-SIGNALS-SURVIVAL",
  "REACTOME-TNFR1-INDUCED-NFKAPPAB-SIGNALING-PATHWAY",
  "REACTOME-REGULATION-OF-TNFR1-SIGNALING",
  "REACTOME-TNFR1-MEDIATED-CERAMIDE-PRODUCTION",
  "REACTOME-IRAK2-MEDIATED-ACTIVATION-OF-TAK1-COMPLEX",
  "REACTOME-IRAK1-RECRUITS-IKK-COMPLEX",
  "REACTOME-TAK1-ACTIVATES-NFKB-BY-PHOSPHORYLATION-AND-ACTIVATION-OF-IKKS-COMPLEX",
  "REACTOME-IKBA-VARIANT-LEADS-TO-EDA-ID",
  "REACTOME-IKK-COMPLEX-RECRUITMENT-MEDIATED-BY-RIP1",
  "REACTOME-TICAM1-RIP1-MEDIATED-IKK-COMPLEX-RECRUITMENT",
  "REACTOME-TICAM1-TRAF6-DEPENDENT-INDUCTION-OF-TAK1-COMPLEX",
  "REACTOME-TRAF6-MEDIATED-INDUCTION-OF-TAK1-COMPLEX-WITHIN-TLR4-COMPLEX",
  "REACTOME-P75NTR-SIGNALS-VIA-NF-KB",
  "REACTOME-P75NTR-RECRUITS-SIGNALLING-COMPLEXES"
 ),
 
 "TGF-beta / EMT / Fibrosis" = c(
  "REACTOME-TGF-BETA-RECEPTOR-SIGNALING-IN-EMT-EPITHELIAL-TO-MESENCHYMAL-TRANSITION",
  "REACTOME-TGF-BETA-RECEPTOR-SIGNALING-ACTIVATES-SMADS",
  "REACTOME-DOWNREGULATION-OF-TGF-BETA-RECEPTOR-SIGNALING",
  "REACTOME-DOWNREGULATION-OF-SMAD2-3-SMAD4-TRANSCRIPTIONAL-ACTIVITY",
  "REACTOME-SMAD2-SMAD3-SMAD4-HETEROTRIMER-REGULATES-TRANSCRIPTION"
 ),
 
 "Innate Immune / IRF / TLR Signaling" = c(
  "REACTOME-TICAM1-DEPENDENT-ACTIVATION-OF-IRF3-IRF7",
  "REACTOME-ACTIVATION-OF-IRF3-IRF7-MEDIATED-BY-TBK1-IKK-EPSILON",
  "REACTOME-TRAF6-MEDIATED-IRF7-ACTIVATION-IN-TLR7-8-OR-9-SIGNALING",
  "REACTOME-REGULATION-OF-INNATE-IMMUNE-RESPONSES-TO-CYTOSOLIC-DNA",
  "REACTOME-NEGATIVE-REGULATORS-OF-DDX58-IFIH1-SIGNALING",
  "REACTOME-REGULATION-OF-BACH1-ACTIVITY",
  "REACTOME-CYTOPROTECTION-BY-HMOX1"
 ),
 
 "MAPK / Stress Kinase Signaling" = c(
  "REACTOME-JNK-C-JUN-KINASES-PHOSPHORYLATION-AND-ACTIVATION-MEDIATED-BY-ACTIVATED-HUMAN-TAK1",
  "REACTOME-ACTIVATED-TAK1-MEDIATES-P38-MAPK-ACTIVATION",
  "REACTOME-MAP3K8-TPL2-DEPENDENT-MAPK1-3-ACTIVATION",
  "REACTOME-NEGATIVE-REGULATION-OF-MAPK-PATHWAY",
  "REACTOME-ACTIVATION-OF-THE-AP-1-FAMILY-OF-TRANSCRIPTION-FACTORS"
 ),
 
 "IL-12 / JAK-STAT Signaling" = c(
  "REACTOME-INTERLEUKIN-12-SIGNALING",
  "REACTOME-INTERLEUKIN-12-FAMILY-SIGNALING",
  "REACTOME-GENE-AND-PROTEIN-EXPRESSION-BY-JAK-STAT-SIGNALING-AFTER-INTERLEUKIN-12-STIMULATION"
 ),
 
 "Antigen Presentation / BCR" = c(
  "REACTOME-ANTIGEN-PRESENTATION-FOLDING-ASSEMBLY-AND-PEPTIDE-LOADING-OF-CLASS-I-MHC",
  "REACTOME-ANTIGEN-PROCESSING-CROSS-PRESENTATION",
  "REACTOME-CROSS-PRESENTATION-OF-PARTICULATE-EXOGENOUS-ANTIGENS-PHAGOSOMES",
  "REACTOME-NEF-MEDIATED-DOWNREGULATION-OF-MHC-CLASS-I-COMPLEX-CELL-SURFACE-EXPRESSION",
  "REACTOME-PD-1-SIGNALING",
  "REACTOME-DOWNSTREAM-SIGNALING-EVENTS-OF-B-CELL-RECEPTOR-BCR"
 ),
 
 "Autophagy / Mitophagy / UPR" = c(
  "REACTOME-CHAPERONE-MEDIATED-AUTOPHAGY",
  "REACTOME-PINK1-PRKN-MEDIATED-MITOPHAGY",
  "REACTOME-MITOPHAGY",
  "REACTOME-PEXOPHAGY",
  "REACTOME-AGGREPHAGY",
  "REACTOME-LATE-ENDOSOMAL-MICROAUTOPHAGY",
  "REACTOME-ATF6-ATF6-ALPHA-ACTIVATES-CHAPERONE-GENES",
  "REACTOME-ATF6-ATF6-ALPHA-ACTIVATES-CHAPERONES",
  "REACTOME-HSF1-ACTIVATION",
  "REACTOME-HSF1-DEPENDENT-TRANSACTIVATION",
  "REACTOME-CALNEXIN-CALRETICULIN-CYCLE",
  "REACTOME-N-GLYCAN-TRIMMING-IN-THE-ER-AND-CALNEXIN-CALRETICULIN-CYCLE",
  "REACTOME-ER-QUALITY-CONTROL-COMPARTMENT-ERQC"
 ),
 
 "Metabolism / Mitochondria" = c(
  "REACTOME-RESPIRATORY-ELECTRON-TRANSPORT",
  "REACTOME-RESPIRATORY-ELECTRON-TRANSPORT-ATP-SYNTHESIS-BY-CHEMIOSMOTIC-COUPLING-AND-HEAT-PRODUCTION-BY-UNCOUPLING-PROTEINS",
  "REACTOME-FORMATION-OF-ATP-BY-CHEMIOSMOTIC-COUPLING",
  "REACTOME-THE-CITRIC-ACID-TCA-CYCLE-AND-RESPIRATORY-ELECTRON-TRANSPORT",
  "REACTOME-COMPLEX-I-BIOGENESIS",
  "REACTOME-CRISTAE-FORMATION",
  "REACTOME-GLYCOGEN-METABOLISM",
  "REACTOME-GLYCOGEN-SYNTHESIS",
  "REACTOME-GLYCOGEN-STORAGE-DISEASES",
  "REACTOME-CELLULAR-RESPONSE-TO-STARVATION",
  "REACTOME-METABOLISM-OF-AMINO-ACIDS-AND-DERIVATIVES",
  "REACTOME-SELENOAMINO-ACID-METABOLISM"
 ),
 
 "Translation / RNA Processing" = c(
  "REACTOME-EUKARYOTIC-TRANSLATION-ELONGATION",
  "REACTOME-EUKARYOTIC-TRANSLATION-INITIATION",
  "REACTOME-TRANSLATION",
  "REACTOME-NONSENSE-MEDIATED-DECAY-NMD",
  "REACTOME-METABOLISM-OF-RNA",
  "REACTOME-RRNA-PROCESSING",
  "REACTOME-AUF1-HNRNP-D0-BINDS-AND-DESTABILIZES-MRNA",
  "REACTOME-REGULATION-OF-MRNA-STABILITY-BY-PROTEINS-THAT-BIND-AU-RICH-ELEMENTS",
  "REACTOME-RESPONSE-OF-EIF2AK4-GCN2-TO-AMINO-ACID-DEFICIENCY",
  "REACTOME-ACTIVATION-OF-THE-MRNA-UPON-BINDING-OF-THE-CAP-BINDING-COMPLEX-AND-EIFS-AND-SUBSEQUENT-BINDING-TO-43S",
  "REACTOME-ATTENUATION-PHASE"
 )
)

# Flatten to a named vector: pathway -> category
pathway_category_map <- stack(pathway_groups)
colnames(pathway_category_map) <- c("pathway", "category")

# Sanity check
cat("Total pathways mapped:", nrow(pathway_category_map), "\n")
cat("Categories:\n")
print(table(pathway_category_map$category))

# Compute mean UCell score per BAFFreceptorStatus group
# Get pathways that are both mapped AND present in ucell_df
available_pathways <- intersect(pathway_category_map$pathway, colnames(ucell_df))
cat("Pathways mapped and present in ucell_df:", length(available_pathways), "\n")

# Subset ucell_df to available pathways + metadata
ucell_sub <- ucell_df[, c(available_pathways, "BAFFreceptorStatus")]

# Compute mean per group
mean_mat <- ucell_sub %>%
 group_by(BAFFreceptorStatus) %>%
 summarise(across(all_of(available_pathways), mean, na.rm = TRUE)) %>%
 column_to_rownames("BAFFreceptorStatus") %>%
 t() %>%
 as.data.frame()

cat("Matrix dims (pathways x groups):", dim(mean_mat), "\n")
print(head(mean_mat))

# Z-score + build row annotation + heatmap
# Z-score across the 2 groups (per pathway)
mean_mat_z <- t(scale(t(mean_mat)))

# Build row annotation from pathway_category_map
row_anno_df <- pathway_category_map %>%
 filter(pathway %in% rownames(mean_mat_z)) %>%
 arrange(match(pathway, rownames(mean_mat_z))) %>%
 column_to_rownames("pathway")

# Reorder mean_mat_z rows to match annotation order
mean_mat_z <- mean_mat_z[rownames(row_anno_df), ]

# Define category colours
n_cats <- length(unique(row_anno_df$category))
cat_colours <- c(RColorBrewer::brewer.pal(8, "Set2"), "#8B4513")
category_cols <- setNames(cat_colours[1:n_cats], unique(row_anno_df$category))

# Clean row names (strip REACTOME- prefix)
rownames(mean_mat_z) <- gsub("^REACTOME-", "", rownames(mean_mat_z))
row_anno_df_plot <- row_anno_df
rownames(row_anno_df_plot) <- gsub("^REACTOME-", "", rownames(row_anno_df_plot))

row_anno <- rowAnnotation(
 Category = row_anno_df_plot$category,
 col = list(Category = category_cols),
 show_annotation_name = FALSE,
 show_legend = TRUE
)

ht <- Heatmap(
 mean_mat_z,
 name = "Z-score",
 column_title = "UCell Pathway Enrichment: HIGH vs LOW expressing BAFF receptors B cells",
 column_title_gp = gpar(fontsize = 14, fontface = "bold"),
 column_title_side = "top",
 col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
 row_split = row_anno_df_plot$category,
 cluster_row_slices = FALSE,
 cluster_rows = FALSE,
 cluster_columns = FALSE,
 left_annotation = row_anno,      # colour bar on LEFT
 show_row_names = TRUE,
 row_names_side = "left",         # names on LEFT before colour bar
 row_title = NULL,
 show_column_names = TRUE,
 column_names_side = "bottom",
 row_names_gp = gpar(fontsize = 12),
 row_names_max_width = unit(12, "cm"),
 border = TRUE,
 width = unit(3, "cm"),
 row_gap = unit(2, "mm"),
 border_gp = gpar(col = "grey60", lwd = 0.5),
 cell_fun = function(j, i, x, y, width, height, fill) {
  grid.rect(x, y, width, height,
            gp = gpar(col = "grey85", fill = NA, lwd = 0.4))
 }
)

png("/data/home/hdx044/plots/BAFF/enrichPathwayUcellHIGH_vs_LOW_BAFFreceptor.png",
    width = 22, height = 20, units = "in", res = 300)
draw(ht,
     merge_legend = TRUE,
     padding = unit(c(2, 2, 2, 2), "mm"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

#### Significant Ucell pathway ####
# Stats: Wilcoxon per pathway, BH correction

pathway_cols <- available_pathways
stat_results <- data.frame(
 pathway = pathway_cols,
 p_value = NA_real_,
 p_adj   = NA_real_
)

for (i in seq_along(pathway_cols)) {
 pw <- pathway_cols[i]
 high <- ucell_df[ucell_df$BAFFreceptorStatus == "High", pw]
 low  <- ucell_df[ucell_df$BAFFreceptorStatus == "Low",  pw]
 stat_results$p_value[i] <- wilcox.test(high, low, exact = FALSE)$p.value
}

stat_results$p_adj <- p.adjust(stat_results$p_value, method = "BH")

# Significance stars
stat_results$stars <- ifelse(stat_results$p_adj < 0.001, "***",
                             ifelse(stat_results$p_adj < 0.01,  "**",
                                    ifelse(stat_results$p_adj < 0.05,  "*", "")))

# Sanity check
cat("Significant pathways (FDR < 0.05):", sum(stat_results$p_adj < 0.05, na.rm = TRUE), "\n")
print(stat_results[stat_results$stars != "", ])

# Add effect size filter: only star if mean difference > threshold
mean_diff <- rowMeans(mean_mat[, "High", drop=FALSE]) - rowMeans(mean_mat[, "Low", drop=FALSE])

stat_results$mean_diff <- mean_diff[stat_results$pathway]
stat_results$stars_filtered <- ifelse(
 stat_results$p_adj < 0.05 & abs(stat_results$mean_diff) > 0.05, 
 stat_results$stars, 
 ""
)
cat("Pathways with FDR<0.05 AND |mean diff|>0.05:", 
    sum(stat_results$stars_filtered != ""), "\n")

# Check distribution of mean differences
summary(abs(stat_results$mean_diff))

stat_results[stat_results$stars_filtered != "", c("pathway", "mean_diff", "p_adj", "stars_filtered")]

#### Highlight pathways with coloured border in cell_fun ####

highlight_paths <- gsub("^REACTOME-", "", c(
 "REACTOME-PEXOPHAGY",
 "REACTOME-ATF6-ATF6-ALPHA-ACTIVATES-CHAPERONE-GENES",
 "REACTOME-ATF6-ATF6-ALPHA-ACTIVATES-CHAPERONES",
 "REACTOME-FORMATION-OF-ATP-BY-CHEMIOSMOTIC-COUPLING",
 "REACTOME-SELENOAMINO-ACID-METABOLISM",
 "REACTOME-EUKARYOTIC-TRANSLATION-ELONGATION",
 "REACTOME-EUKARYOTIC-TRANSLATION-INITIATION",
 "REACTOME-NONSENSE-MEDIATED-DECAY-NMD",
 "REACTOME-RESPONSE-OF-EIF2AK4-GCN2-TO-AMINO-ACID-DEFICIENCY",
 "REACTOME-ACTIVATION-OF-THE-MRNA-UPON-BINDING-OF-THE-CAP-BINDING-COMPLEX-AND-EIFS-AND-SUBSEQUENT-BINDING-TO-43S"
))

highlight_rows <- which(rownames(mean_mat_z) %in% highlight_paths)

ht <- Heatmap(
 mean_mat_z,
 name = "Z-score",
 col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
 row_split = row_anno_df_plot$category,
 cluster_row_slices = FALSE,
 cluster_rows = FALSE,
 cluster_columns = FALSE,
 left_annotation = row_anno,
 show_row_names = TRUE,
 row_names_side = "left",
 row_title = NULL,
 show_column_names = TRUE,
 column_names_side = "bottom",
 row_names_gp = gpar(fontsize = 12),
 row_names_max_width = unit(12, "cm"),
 border = TRUE,
 width = unit(3, "cm"),
 row_gap = unit(2, "mm"),
 border_gp = gpar(col = "grey60", lwd = 0.5),
 cell_fun = function(j, i, x, y, width, height, fill) {
  grid.rect(x, y, width, height,
            gp = gpar(col = "grey85", fill = NA, lwd = 0.4))
  if (i %in% highlight_rows) {
   grid.rect(x, y, width, height,
             gp = gpar(col = "blue", fill = NA, lwd = 2))
  }
 }
)

png("/data/home/hdx044/plots/BAFF/enrichPathwayUcellHIGH_vs_LOW_BAFFreceptor.png",
    width = 22, height = 20, units = "in", res = 300)
draw(ht,
     merge_legend = TRUE,
     padding = unit(c(15, 2, 2, 2), "mm"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
grid.text("UCell Pathway Enrichment: BAFF-Receptor High vs Low B cells",
          x = unit(2, "mm"), y = unit(1, "npc") - unit(5, "mm"),
          just = "left",
          gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()


#### Create cellchat obj for baff receptor ####

# Split Seurat Object Based on BAFF Receptor Expression Levels
# Creates two objects: HIGH receptor (>0.25) vs LOW receptor (≤0.25) in specific clusters
# All other clusters remain unchanged in both objects

# BAFF receptor genes (human gene names)

baff_receptors <- c("TNFRSF13C", "TNFRSF13B", "TNFRSF17")  # BAFFR, TACI, BCMA

# Clusters where we want to filter by receptor expression
target_clusters <- c(8, 19, 22, 27, 28, 29)

# Expression threshold
expression_threshold <- 0.25

# Cluster column name in metadata (adjust if different)
cluster_column <- "seurat_clusters"  # or "cluster", "celltype", etc.

cat("Seurat object loaded successfully!\n")
cat(paste0("Total cells: ", ncol(LIVER), "\n"))
cat(paste0("Total genes: ", nrow(LIVER), "\n\n"))

# Check if cluster column exists
if (!cluster_column %in% colnames(LIVER@meta.data)) {
 cat("Available metadata columns:\n")
 print(colnames(LIVER@meta.data))
 stop(paste0("Error: '", cluster_column, "' not found in metadata!\n",
             "Please update the 'cluster_column' parameter."))
}

# CHECK BAFF RECEPTOR GENES
cat("Checking BAFF receptor genes...\n")
receptors_present <- baff_receptors[baff_receptors %in% rownames(LIVER)]
receptors_missing <- baff_receptors[!baff_receptors %in% rownames(LIVER)]

if (length(receptors_present) == 0) {
 stop("Error: None of the BAFF receptor genes found in the object!\n",
      "Available genes: ", paste(head(rownames(LIVER), 20), collapse=", "))
}

cat(paste0("Found receptors: ", paste(receptors_present, collapse=", "), "\n"))
if (length(receptors_missing) > 0) {
 cat(paste0("Missing receptors: ", paste(receptors_missing, collapse=", "), "\n"))
}
cat("\n")

# EXTRACT BAFF RECEPTOR EXPRESSION
cat("Extracting BAFF receptor expression...\n")

# Get expression matrix for all receptor genes present
receptor_expression <- FetchData(LIVER, vars = receptors_present)

# Calculate maximum expression across all BAFF receptors for each cell
# This identifies cells expressing ANY receptor above threshold
receptor_expression$max_receptor <- apply(receptor_expression[, receptors_present, drop=FALSE], 
                                          1, max)

# Add cluster information
receptor_expression$cluster <- LIVER@meta.data[[cluster_column]]

# Summary statistics
cat("\nReceptor expression summary (max across all receptors):\n")
print(summary(receptor_expression$max_receptor))

# IDENTIFY CELLS FOR EACH OBJECT

# Get all cell barcodes
all_cells <- colnames(LIVER)

# Initialize vectors to store selected cells
cells_high_receptor <- c()
cells_low_receptor <- c()

# Process each cluster
for (clust in unique(LIVER@meta.data[[cluster_column]])) {
 
 # Get cells in this cluster
 cluster_cells <- rownames(LIVER@meta.data[LIVER@meta.data[[cluster_column]] == clust, ])
 
 if (clust %in% target_clusters) {
  # TARGET CLUSTERS: Filter by receptor expression
  
  # Get receptor expression for this cluster
  cluster_receptor_expr <- receptor_expression[cluster_cells, ]
  
  # HIGH receptor cells (> 0.25)
  high_cells <- rownames(cluster_receptor_expr[cluster_receptor_expr$max_receptor > expression_threshold, ])
  
  # LOW receptor cells (≤ 0.25)
  low_cells <- rownames(cluster_receptor_expr[cluster_receptor_expr$max_receptor <= expression_threshold, ])
  
  cells_high_receptor <- c(cells_high_receptor, high_cells)
  cells_low_receptor <- c(cells_low_receptor, low_cells)
  
  cat(sprintf("Cluster %s (TARGET):\n", clust))
  cat(sprintf("  Total cells: %d\n", length(cluster_cells)))
  cat(sprintf("  HIGH receptor (>%.2f): %d cells (%.1f%%)\n", 
              expression_threshold, length(high_cells), 
              length(high_cells)/length(cluster_cells)*100))
  cat(sprintf("  LOW receptor (≤%.2f): %d cells (%.1f%%)\n\n", 
              expression_threshold, length(low_cells), 
              length(low_cells)/length(cluster_cells)*100))
  
 } else {
  # NON-TARGET CLUSTERS: Keep all cells in both objects
  
  cells_high_receptor <- c(cells_high_receptor, cluster_cells)
  cells_low_receptor <- c(cells_low_receptor, cluster_cells)
  
  cat(sprintf("Cluster %s (keep all): %d cells\n", clust, length(cluster_cells)))
 }
}

cat("FINAL CELL COUNTS\n")
cat(sprintf("Original object: %d cells\n", ncol(LIVER)))
cat(sprintf("HIGH receptor object: %d cells\n", length(cells_high_receptor)))
cat(sprintf("LOW receptor object: %d cells\n", length(cells_low_receptor)))
cat("\n")

# CREATE TWO SEURAT OBJECTS

cat("Creating HIGH receptor Seurat object...\n")
seurat_high_receptor <- subset(LIVER, cells = cells_high_receptor)

cat("Creating LOW receptor Seurat object...\n")
seurat_low_receptor <- subset(LIVER, cells = cells_low_receptor)

cat("\nHIGH receptor object summary:\n")
cat(sprintf("  Cells: %d\n", ncol(seurat_high_receptor)))
cat(sprintf("  Genes: %d\n", nrow(seurat_high_receptor)))

cat("\nLOW receptor object summary:\n")
cat(sprintf("  Cells: %d\n", ncol(seurat_low_receptor)))
cat(sprintf("  Genes: %d\n", nrow(seurat_low_receptor)))

# CLUSTER COMPOSITION COMPARISON

# Get cluster compositions
original_comp <- table(LIVER@meta.data[[cluster_column]])
high_comp <- table(seurat_high_receptor@meta.data[[cluster_column]])
low_comp <- table(seurat_low_receptor@meta.data[[cluster_column]])

# Create comparison data frame
comp_df <- data.frame(
 Cluster = names(original_comp),
 Original = as.numeric(original_comp),
 HIGH_receptor = as.numeric(high_comp[names(original_comp)]),
 LOW_receptor = as.numeric(low_comp[names(original_comp)])
)
comp_df[is.na(comp_df)] <- 0

# Add percentages
comp_df$HIGH_pct <- round(comp_df$HIGH_receptor / sum(comp_df$HIGH_receptor) * 100, 2)
comp_df$LOW_pct <- round(comp_df$LOW_receptor / sum(comp_df$LOW_receptor) * 100, 2)

# Mark target clusters
comp_df$Is_Target <- ifelse(as.numeric(as.character(comp_df$Cluster)) %in% target_clusters, 
                            "YES", "NO")

# Reorder columns
comp_df <- comp_df[, c("Cluster", "Is_Target", "Original", "HIGH_receptor", 
                       "HIGH_pct", "LOW_receptor", "LOW_pct")]

cat("Cluster composition:\n")
print(comp_df)

# RECEPTOR EXPRESSION ANALYSIS IN TARGET CLUSTERS
# Analyze each receptor individually in target clusters
for (receptor in receptors_present) {
 cat(paste0("\n", receptor, ":\n"))
 cat(paste0(rep("-", nchar(receptor) + 1), collapse=""), "\n")
 
 receptor_expr <- FetchData(LIVER, vars = receptor)
 receptor_expr$cluster <- LIVER@meta.data[[cluster_column]]
 
 for (clust in target_clusters) {
  cluster_expr <- receptor_expr[receptor_expr$cluster == clust, receptor]
  n_expressing <- sum(cluster_expr > expression_threshold)
  pct_expressing <- round(n_expressing / length(cluster_expr) * 100, 2)
  mean_expr <- round(mean(cluster_expr), 3)
  
  cat(sprintf("  Cluster %s: %.1f%% cells >%.2f (mean=%.3f)\n", 
              clust, pct_expressing, expression_threshold, mean_expr))
 }
}

# SAVE OBJECTS

# Save Seurat objects
high_file <- file.path('/data/Blizard-AlazawiLab/rk/cellchat/BAFF', "liver_seurat_HIGH_receptor.rds")
low_file <- file.path('/data/Blizard-AlazawiLab/rk/cellchat/BAFF', "liver_seurat_LOW_receptor.rds")

cat(paste0("Saving HIGH receptor object to: ", high_file, "\n"))
saveRDS(seurat_high_receptor, file = high_file)

cat(paste0("Saving LOW receptor object to: ", low_file, "\n"))
saveRDS(seurat_low_receptor, file = low_file)

# Save composition table
comp_file <- file.path('/data/Blizard-AlazawiLab/rk/cellchat/BAFF', "cluster_composition.csv")
write.csv(comp_df, comp_file, row.names = FALSE)
cat(paste0("Saved composition table to: ", comp_file, "\n"))


#### Load FNA data ####
SeuObj <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_LIVER_processed.rds")

# BAFF and its receptor genes (human gene names)
baff_receptors <- c("TNFSF13B", "TNFRSF13C", "TNFRSF13B", "TNFRSF17")  # BAFF, BAFFR, TACI, BCMA

DotPlot2(SeuObj, features = baff_receptors, split.by = "seurat_clusters")

fnaBAFF <- subset(SeuObj, subset = seurat_clusters %in% c("6", "9", "12", "21"))

fnaBAFF$PatientID <- sub(".* ", "", fnaBAFF$VisitLabel)
table(fnaBAFF$PatientID)

# Get all unique Visit Labels
visits <- unique(fnaBAFF$VisitLabel)

# Extract patient IDs
patient_ids <- sub(".* ", "", visits)

# Extract visit type (baseline / followup)
visit_type <- ifelse(grepl("^baseline", visits), "baseline", "followup")

df <- data.frame(
 VisitLabel = visits,
 PatientID = patient_ids,
 VisitType = visit_type,
 stringsAsFactors = FALSE
)

# Order: baseline BEFORE followup, within each patient
df <- df[order(df$PatientID, df$VisitType), ]

# Apply ordered factor to your object
fnaBAFF$VisitLabel <- factor(fnaBAFF$VisitLabel, levels = df$VisitLabel)
levels(fnaBAFF$VisitLabel)


ggp <- DotPlot2(
 fnaBAFF,
 features = baff_receptors,
 group.by = "VisitLabel"
)

ggp

ggsave("/data/home/hdx044/plots/screpertoire/liver/FNA/BCR/FNAdiversityBCR.png",
       ggp, width = 5, height = 5, dpi = 300)




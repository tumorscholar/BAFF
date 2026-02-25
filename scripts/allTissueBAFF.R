# Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(SeuratExtend)

# Load seurat object
SeuObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/SeuObjx.rds')

DefaultAssay(SeuObj) <- "RNA"

baff_signaling_genes <- c(
 "TNFRSF13C",  # BAFF-R (BAFF receptor)
 "TNFRSF13B",  # TACI (Transmembrane activator and CAML interactor)
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

# BAFF and its receptors clusters
liverBaff <- subset(LIVER, subset = cluster %in% c("C3", "C5", "C8", "C10", "C11", "C13", "C14", "C19", "C22", "C23", "C27", "C29"))

saveRDS(liverBaff, '/data/Blizard-AlazawiLab/rk/seurat/liverBaff.rds')

#### load liver BAFF obj ####
liverBaff <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/liverBaff.rds')

Idents(liverBaff) <- "Stage"
DefaultAssay(liverBaff) <- "RNA"

# Arrange Stage in the desired order
liverBaff$Stage <- factor(
 liverBaff$Stage,
 levels = c("Healthy", "No_fibrosis", "Fibrosis")
)

ggp <- DotPlot2(liverBaff, features = baff_signaling_genes, color_scheme = "BuRd", group.by = "Stage") +
 ggtitle("BAFF and its receptors") +
 guides(
  color = guide_colorbar(order = 1),  # Average Expression on top
  size = guide_legend(order = 2))      # Percent Expressed below

ggp 

# Save high-resolution SVG
ggsave(
 filename = "LiverBAFFandReceptorsStage.svg",
 plot = ggp,
 device = "svg",
 path = "/data/home/hdx044/plots/BAFF",
 width = 3,     
 height = 4,
 units = "in"
)


baff_signaling_ADTs <- c(
 "Hu.CD268",  #  TNFRSF13C also BAFF-R (BAFF receptor)
 "Hu.CD267"  #  TNFRSF17 also BCMA (B cell maturation antigen)
)

DefaultAssay(liverBaff) <- "ADTonly"

ggp <- DotPlot2(liverBaff, features = baff_signaling_ADTs, color_scheme = "BuRd", group.by = "Stage") +
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

#### Create cellchat obj ####
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

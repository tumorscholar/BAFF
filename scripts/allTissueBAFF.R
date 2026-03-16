# Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(SeuratExtend)
library(escape)

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

#### Pathway analysis ####
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

#### load liver BAFF obj ####
liverBaffreceptorObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/liverBaffreceptorObj.rds')

library(dplyr)
library(tidyr)
library(pheatmap)

ucell_mat <- GetAssayData(
 liverBaffreceptorObj,
 assay = "escape.UCell",
 slot = "data"
)

ucell_df <- as.data.frame(
 t(as.matrix(ucell_mat)),
 check.names = FALSE
)

# USE CELL TYPE NAMES HERE
ucell_df$cluster <- factor(liverBaffreceptorObj$cell_type_with_cluster)

pathway_cluster_means <- ucell_df %>%
 group_by(.data$cluster) %>%
 summarise(across(everything(), ~mean(.x)), .groups = "drop") %>%
 relocate(cluster)

long_means <- pathway_cluster_means %>%
 pivot_longer(
  -cluster,
  names_to = "pathway",
  values_to = "mean_score"
 )

write.csv(
 long_means,
 file = "/data/home/hdx044/files/BAFF/All_pathways_by_celltype_BAFF.csv",
 row.names = FALSE
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

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

# EXPANDED PATHWAY FILTERING

# Define comprehensive pathway categories relevant to MASH and B cell biology

pathway_categories <- list(
 
 # FIBROSIS-RELATED
 fibrosis = c(
  "TGF", "TGFB", "PDGF", "WNT", "BMP", "VEGF",
  "COLLAGEN", "ECM", "EXTRACELLULAR.MATRIX", "MATRIX",
  "FIBRONECTIN", "EMT", "EPITHELIAL", "FIBROSIS", "FIBROTIC"
 ),
 
 # B CELL ACTIVATION & SIGNALING
 activation = c(
  "B.CELL", "BCR", "B.LYMPHOCYTE",
  "CD40", "BAFF", "TNFSF13", "APRIL",
  "ACTIVATION", "ANTIGEN.RECEPTOR",
  "NF.KAPPA.B", "NFKB", "PI3K", "AKT",
  "MTOR", "JAK", "STAT"
 ),
 
 # PROLIFERATION & CELL CYCLE
 proliferation = c(
  "PROLIFERATION", "CELL.CYCLE", "MITOSIS", "MITOTIC",
  "G1.S", "G2.M", "CYCLIN", "CDK",
  "DNA.REPLICATION", "S.PHASE", "M.PHASE",
  "E2F", "RB1", "P53"
 ),
 
 # INFLAMMATION & CYTOKINES
 inflammation = c(
  "INFLAMMATION", "INFLAMMATORY",
  "TNF", "TNFSF", "IL.1", "IL.6", "IL.10",
  "INTERFERON", "IFN", "CYTOKINE",
  "CHEMOKINE", "CCL", "CXCL"
 ),
 
 # IMMUNE RESPONSE
 immune = c(
  "IMMUNE", "IMMUNOLOGICAL",
  "COMPLEMENT", "C3", "C5",
  "MHC", "ANTIGEN.PRESENTATION",
  "COSTIMULATORY", "CD80", "CD86"
 ),
 
 # APOPTOSIS & SURVIVAL
 survival = c(
  "APOPTOSIS", "APOPTOTIC", "CELL.DEATH",
  "SURVIVAL", "BCL", "CASPASE",
  "INTRINSIC.APOPTOSIS", "EXTRINSIC.APOPTOSIS"
 ),
 
 # METABOLISM
 metabolism = c(
  "METABOLISM", "METABOLIC", "GLYCOLYSIS",
  "OXIDATIVE.PHOSPHORYLATION", "OXPHOS",
  "FATTY.ACID", "LIPID", "GLUCOSE"
 ),
 
 # PLASMA CELL DIFFERENTIATION
 plasma_cell = c(
  "PLASMA.CELL", "ANTIBODY", "IMMUNOGLOBULIN",
  "PRDM1", "BLIMP", "XBP1", "IRF4",
  "SECRETION", "ER.STRESS", "UPR"
 )
)

# Create combined pattern for each category
category_patterns <- lapply(pathway_categories, function(keywords) {
 paste(keywords, collapse = "|")
})

# FILTER PATHWAYS BY CATEGORY

cat("\n========== PATHWAY FILTERING BY CATEGORY ==========\n\n")

pathway_results <- list()

for (cat_name in names(category_patterns)) {
 
 pattern <- category_patterns[[cat_name]]
 
 pathways_cat <- pathway_max_scores %>%
  filter(
   max_score > UCELL_MODERATE,  # Still use 0.15 cutoff
   grepl(pattern, pathway, ignore.case = TRUE)
  ) %>%
  arrange(desc(max_score))
 
 pathway_results[[cat_name]] <- pathways_cat
 
 cat(sprintf("%s pathways (>0.15): %d\n", 
             toupper(cat_name), nrow(pathways_cat)))
 
 if (nrow(pathways_cat) > 0) {
  cat("  Top 5:\n")
  for (i in 1:min(5, nrow(pathways_cat))) {
   cat(sprintf("    %.3f - %s\n", 
               pathways_cat$max_score[i],
               substr(pathways_cat$pathway[i], 1, 70)))
  }
 }
 cat("\n")
}

# COMBINE ALL RELEVANT PATHWAYS

# Combine all categories (removing duplicates)
pathways_combined <- unique(c(
 pathway_results$fibrosis$pathway,
 pathway_results$activation$pathway,
 pathway_results$proliferation$pathway,
 pathway_results$inflammation$pathway,
 pathway_results$immune$pathway,
 pathway_results$survival$pathway,
 pathway_results$metabolism$pathway,
 pathway_results$plasma_cell$pathway
))

cat(sprintf("TOTAL UNIQUE PATHWAYS (all categories): %d\n\n", 
            length(pathways_combined)))

# ANNOTATE PATHWAYS WITH CATEGORIES

pathway_annotations <- pathway_max_scores %>%
 filter(pathway %in% pathways_combined) %>%
 mutate(
  category = case_when(
   grepl(category_patterns$fibrosis, pathway, ignore.case = TRUE) ~ "Fibrosis",
   grepl(category_patterns$activation, pathway, ignore.case = TRUE) ~ "B cell activation",
   grepl(category_patterns$proliferation, pathway, ignore.case = TRUE) ~ "Proliferation",
   grepl(category_patterns$inflammation, pathway, ignore.case = TRUE) ~ "Inflammation",
   grepl(category_patterns$plasma_cell, pathway, ignore.case = TRUE) ~ "Plasma cell",
   grepl(category_patterns$immune, pathway, ignore.case = TRUE) ~ "Immune response",
   grepl(category_patterns$survival, pathway, ignore.case = TRUE) ~ "Survival/Apoptosis",
   grepl(category_patterns$metabolism, pathway, ignore.case = TRUE) ~ "Metabolism",
   TRUE ~ "Other"
  )
 ) %>%
 arrange(category, desc(max_score))

# CREATE COMPREHENSIVE HEATMAP

heatmap_data_full <- pathway_cluster_means %>%
 select(cluster, all_of(pathways_combined)) %>%
 as.data.frame()

rownames(heatmap_data_full) <- heatmap_data_full$cluster
heatmap_matrix_full <- t(as.matrix(heatmap_data_full[, -1]))

# Create category annotation for heatmap
pathway_category_anno <- pathway_annotations %>%
 select(pathway, category) %>%
 arrange(match(pathway, rownames(heatmap_matrix_full)))

library(dplyr)
library(stringr)
library(pheatmap)
library(viridis)

# make sure matrix rownames match "pathway"
pathway_category_anno <- pathway_annotations %>%
 select(pathway, category)

categories_to_plot <- c(
 "Fibrosis",
 "B cell activation",
 "Proliferation",
 "Inflammation",
 "Plasma cell",
 "Immune response",
 "Survival/Apoptosis",
 "Metabolism"
)

# Directory to save heatmaps
outdir <- "/data/home/hdx044/plots/BAFF/category_heatmaps/"
dir.create(outdir, showWarnings = FALSE)

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(grid)
library(officer)
library(rvg)

outdir <- "/data/home/hdx044/plots/BAFF/category_heatmaps/"
dir.create(outdir, showWarnings = FALSE)

# Single PPT file path
ppt_all <- file.path(outdir, "all_categories_heatmaps.pptx")
doc <- read_pptx()  # <-- open once

for (cat in categories_to_plot) {
 # Safe category label for messages only
 safe_cat <- gsub("[^A-Za-z0-9]+", "_", cat)
 
 # Subset matrix for this category
 cat_paths <- pathway_category_anno %>%
  filter(category == cat) %>%
  pull(pathway)
 
 sub_mat <- heatmap_matrix_full[
  rownames(heatmap_matrix_full) %in% cat_paths, , drop = FALSE
 ]
 if (nrow(sub_mat) == 0) {
  message("Skip empty category: ", cat)
  next
 }
 
 # Color function (fix breaks = colors)
 col_fun <- colorRamp2(
  seq(min(sub_mat), max(sub_mat), length.out = 5),
  viridis::viridis(5)
 )
 
 # Build ComplexHeatmap
 p <- Heatmap(
  sub_mat,
  name = "UCell",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  # Keep full Reactome names; adjust font if needed
  row_names_gp = gpar(fontsize = 6),
  row_names_max_width = unit(9, "cm"),
  
  # Column labels
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 10),
  
  # Control block size (cell height/width equivalents)
  height = unit(nrow(sub_mat) * 3.5, "mm"),
  width  = unit(ncol(sub_mat) * 15,  "mm"),
  
  column_title = paste("Category:", cat)
 )
 
 # ---- Add a slide for this category and put the heatmap on it
 doc <- add_slide(doc, layout = "Blank", master = "Office Theme")
 doc <- ph_with(
  doc,
  value = dml(code = { draw(p) }),
  location = ph_location(left = 0.5, top = 0.5, width = 9, height = 5)
 )
 
 message("Added to PPT (slide): ", safe_cat)
}

# Save ONE PPT containing all slides
print(doc, target = ppt_all)
message("Saved single PPT: ", ppt_all)


write.csv(
 heatmap_matrix_full,
 file = "/data/home/hdx044/files/BAFF/pathway_cluster_heatmap_metadata.csv",
 row.names = TRUE
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

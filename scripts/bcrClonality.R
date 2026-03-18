# ============================================================================
# BCR Clonality and BAFF Receptor Analysis
# ============================================================================
# This script analyzes:
# 1. BCR V(D)J gene usage patterns
# 2. Somatic hypermutation markers
# 3. Correlation between clonality, receptor status, and effector functions

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# ============================================================================
# SECTION 1: BCR REPERTOIRE ANALYSIS
# ============================================================================

cat("=" %R% 80, "\n")
cat("SECTION 1: BCR V(D)J GENE USAGE ANALYSIS\n")
cat("=" %R% 80, "\n\n")

# Assuming your Seurat object has BCR data in metadata
# Common column names: IGH_v_gene, IGH_d_gene, IGH_j_gene, IGH_c_gene
# Or: CTaa, CTnt (for clonotype)

# Check what BCR-related columns exist
bcr_columns <- grep("^IG[HKL]|^CTaa|^CTnt|clone", 
                    colnames(seuratObj@meta.data), 
                    value = TRUE, ignore.case = TRUE)

cat("Available BCR-related columns:\n")
print(bcr_columns)
cat("\n")

# ============================================================================
# 1.1: V(D)J Gene Usage in Hyperexpanded Clones
# ============================================================================

cat("1.1: Analyzing V(D)J gene usage in hyperexpanded clones\n")
cat(rep("-", 80), "\n\n")

# Filter to cells with clonotype information and receptor status
bcr_data <- seuratObj@meta.data %>%
 filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize))

# Identify hyperexpanded clones in LOW receptor cells
hyperexpanded_low <- bcr_data %>%
 filter(BAFFreceptorStatus == "Low", 
        cloneSize == "Hyperexpanded (100 < X <= 500)")

cat(sprintf("Number of cells in hyperexpanded LOW receptor clones: %d\n", 
            nrow(hyperexpanded_low)))

# If you have V gene information (adjust column name as needed)
# Example with IGH_v_gene or similar
if ("IGH_v_gene" %in% colnames(bcr_data)) {
 
 # V gene usage in hyperexpanded LOW vs HIGH
 v_gene_comparison <- bcr_data %>%
  filter(cloneSize %in% c("Hyperexpanded (100 < X <= 500)", 
                          "Large (20 < X <= 100)")) %>%
  group_by(BAFFreceptorStatus, IGH_v_gene) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(BAFFreceptorStatus) %>%
  mutate(freq = n / sum(n) * 100) %>%
  arrange(BAFFreceptorStatus, desc(freq))
 
 # Top 10 V genes per group
 cat("\nTop 10 V genes in expanded clones:\n\n")
 cat("LOW receptor:\n")
 print(v_gene_comparison %>% filter(BAFFreceptorStatus == "Low") %>% head(10))
 
 cat("\nHIGH receptor:\n")
 print(v_gene_comparison %>% filter(BAFFreceptorStatus == "High") %>% head(10))
 
 # Plot V gene usage
 top_v_genes <- v_gene_comparison %>%
  group_by(IGH_v_gene) %>%
  summarise(total = sum(n)) %>%
  arrange(desc(total)) %>%
  head(20) %>%
  pull(IGH_v_gene)
 
 p_vgene <- v_gene_comparison %>%
  filter(IGH_v_gene %in% top_v_genes) %>%
  ggplot(aes(x = reorder(IGH_v_gene, -freq), y = freq, fill = BAFFreceptorStatus)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "V Gene Usage in Expanded B Cell Clones",
       subtitle = "Hyperexpanded + Large clones",
       x = "V Gene", y = "Frequency (%)",
       fill = "BAFF Receptor") +
  scale_fill_manual(values = c("Low" = "#E41A1C", "High" = "#377EB8"))
 
 ggsave("v_gene_usage_comparison.pdf", p_vgene, width = 12, height = 6)
 cat("\nSaved: v_gene_usage_comparison.pdf\n")
}

# ============================================================================
# 1.2: Shared Clones Analysis
# ============================================================================

cat("\n1.2: Analyzing shared clones across samples/patients\n")
cat(rep("-", 80), "\n\n")

# Assuming you have a sample/patient ID column
if ("orig.ident" %in% colnames(bcr_data) | "patient" %in% colnames(bcr_data)) {
 
 patient_col <- if("patient" %in% colnames(bcr_data)) "patient" else "orig.ident"
 
 # Identify clones present in multiple patients
 if ("CTaa" %in% colnames(bcr_data)) {
  
  clone_sharing <- bcr_data %>%
   filter(!is.na(CTaa), CTaa != "") %>%
   group_by(CTaa) %>%
   summarise(
    n_cells = n(),
    n_patients = n_distinct(!!sym(patient_col)),
    patients = paste(unique(!!sym(patient_col)), collapse = ","),
    receptor_status = paste(unique(BAFFreceptorStatus), collapse = ","),
    .groups = "drop"
   ) %>%
   arrange(desc(n_patients), desc(n_cells))
  
  # Public clones (shared across patients)
  public_clones <- clone_sharing %>%
   filter(n_patients > 1)
  
  cat(sprintf("\nTotal unique clonotypes: %d\n", nrow(clone_sharing)))
  cat(sprintf("Public clones (shared across patients): %d\n", nrow(public_clones)))
  
  if (nrow(public_clones) > 0) {
   cat("\nTop 20 public clones:\n")
   print(public_clones %>% head(20))
   
   write.csv(public_clones, "public_clones_shared.csv", row.names = FALSE)
   cat("\nSaved: public_clones_shared.csv\n")
   
   # Are public clones enriched in LOW or HIGH receptor?
   public_clone_receptor <- bcr_data %>%
    filter(CTaa %in% public_clones$CTaa) %>%
    group_by(BAFFreceptorStatus) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(freq = n / sum(n) * 100)
   
   cat("\nReceptor status of cells in public clones:\n")
   print(public_clone_receptor)
  }
 }
}

# ============================================================================
# SECTION 2: SOMATIC HYPERMUTATION MARKERS
# ============================================================================

cat("\n\n")
cat("=" %R% 80, "\n")
cat("SECTION 2: SOMATIC HYPERMUTATION AND ANTIGEN EXPERIENCE\n")
cat("=" %R% 80, "\n\n")

# ============================================================================
# 2.1: AICDA Expression Analysis
# ============================================================================

cat("2.1: AICDA (activation-induced deaminase) expression\n")
cat(rep("-", 80), "\n\n")

# Check if AICDA is in the dataset
if ("AICDA" %in% rownames(seuratObj)) {
 
 # Extract AICDA expression
 aicda_expr <- FetchData(seuratObj, 
                         vars = c("AICDA", "BAFFreceptorStatus", "cloneSize"))
 
 # Compare by receptor status
 aicda_summary <- aicda_expr %>%
  filter(!is.na(BAFFreceptorStatus)) %>%
  group_by(BAFFreceptorStatus) %>%
  summarise(
   mean_AICDA = mean(AICDA),
   median_AICDA = median(AICDA),
   pct_expressing = sum(AICDA > 0) / n() * 100,
   .groups = "drop"
  )
 
 cat("\nAICDA expression by receptor status:\n")
 print(aicda_summary)
 
 # Statistical test
 aicda_test <- wilcox.test(AICDA ~ BAFFreceptorStatus, 
                           data = aicda_expr %>% filter(!is.na(BAFFreceptorStatus)))
 cat(sprintf("\nWilcoxon test p-value: %.2e\n", aicda_test$p.value))
 
 # By clonality
 aicda_clonality <- aicda_expr %>%
  filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize)) %>%
  group_by(BAFFreceptorStatus, cloneSize) %>%
  summarise(
   mean_AICDA = mean(AICDA),
   pct_expressing = sum(AICDA > 0) / n() * 100,
   n_cells = n(),
   .groups = "drop"
  )
 
 cat("\nAICDA expression by clonality:\n")
 print(aicda_clonality)
 
 # Plot
 p_aicda <- ggplot(aicda_expr %>% filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize)),
                   aes(x = cloneSize, y = AICDA, fill = BAFFreceptorStatus)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "AICDA Expression by Clonality and BAFF Receptor Status",
       x = "Clone Size", y = "AICDA Expression",
       fill = "BAFF Receptor") +
  scale_fill_manual(values = c("Low" = "#E41A1C", "High" = "#377EB8"))
 
 ggsave("aicda_expression_clonality.pdf", p_aicda, width = 12, height = 6)
 cat("\nSaved: aicda_expression_clonality.pdf\n")
 
} else {
 cat("AICDA not found in dataset\n")
}

# ============================================================================
# 2.2: Other Markers of Antigen Experience
# ============================================================================

cat("\n2.2: Additional markers of antigen experience\n")
cat(rep("-", 80), "\n\n")

# Markers to check
antigen_markers <- c(
 "AICDA",      # Somatic hypermutation
 "CD27",       # Memory B cell
 "CD38",       # Activated/plasma cell
 "IRF4",       # Plasma cell differentiation
 "PRDM1",      # BLIMP1 - plasma cell
 "XBP1",       # ER stress/plasma cell
 "MZB1",       # Marginal zone/memory
 "TNFRSF13B",  # TACI
 "TNFRSF17"    # BCMA
)

# Check which are available
available_markers <- antigen_markers[antigen_markers %in% rownames(seuratObj)]
cat("Available antigen experience markers:\n")
print(available_markers)
cat("\n")

if (length(available_markers) > 0) {
 
 # Extract expression
 marker_expr <- FetchData(seuratObj, 
                          vars = c(available_markers, "BAFFreceptorStatus", "cloneSize"))
 
 # Compare by clonality within receptor groups
 marker_comparison <- marker_expr %>%
  filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize)) %>%
  group_by(BAFFreceptorStatus, cloneSize) %>%
  summarise(across(all_of(available_markers), 
                   list(mean = mean, pct = ~sum(. > 0) / n() * 100)),
            n = n(),
            .groups = "drop")
 
 write.csv(marker_comparison, "antigen_experience_markers_by_clonality.csv", row.names = FALSE)
 cat("Saved: antigen_experience_markers_by_clonality.csv\n")
 
 # Focus on hyperexpanded clones
 hyperexpanded_markers <- marker_expr %>%
  filter(!is.na(BAFFreceptorStatus),
         cloneSize == "Hyperexpanded (100 < X <= 500)") %>%
  group_by(BAFFreceptorStatus) %>%
  summarise(across(all_of(available_markers), 
                   list(mean = mean, pct = ~sum(. > 0) / n() * 100)),
            n = n(),
            .groups = "drop")
 
 cat("\nMarkers in HYPEREXPANDED clones:\n")
 print(hyperexpanded_markers)
}

# ============================================================================
# SECTION 3: EFFECTOR FUNCTION vs CLONALITY
# ============================================================================

cat("\n\n")
cat("=" %R% 80, "\n")
cat("SECTION 3: EFFECTOR MOLECULES vs CLONALITY\n")
cat("=" %R% 80, "\n\n")

# ============================================================================
# 3.1: TGFB1 Expression by Clonality
# ============================================================================

cat("3.1: TGFB1 expression in hyperexpanded LOW receptor cells\n")
cat(rep("-", 80), "\n\n")

if ("TGFB1" %in% rownames(seuratObj)) {
 
 tgfb1_data <- FetchData(seuratObj, 
                         vars = c("TGFB1", "BAFFreceptorStatus", "cloneSize"))
 
 # Summary by group
 tgfb1_summary <- tgfb1_data %>%
  filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize)) %>%
  group_by(BAFFreceptorStatus, cloneSize) %>%
  summarise(
   mean_TGFB1 = mean(TGFB1),
   median_TGFB1 = median(TGFB1),
   pct_expressing = sum(TGFB1 > 0) / n() * 100,
   n_cells = n(),
   .groups = "drop"
  )
 
 cat("TGFB1 expression by clonality:\n")
 print(tgfb1_summary)
 
 # Test: Hyperexpanded LOW vs other LOW
 low_hyperexp <- tgfb1_data %>%
  filter(BAFFreceptorStatus == "Low",
         cloneSize == "Hyperexpanded (100 < X <= 500)") %>%
  pull(TGFB1)
 
 low_other <- tgfb1_data %>%
  filter(BAFFreceptorStatus == "Low",
         cloneSize != "Hyperexpanded (100 < X <= 500)") %>%
  pull(TGFB1)
 
 if (length(low_hyperexp) > 0 & length(low_other) > 0) {
  tgfb1_test <- wilcox.test(low_hyperexp, low_other)
  cat(sprintf("\nHyperexpanded LOW vs other LOW: p = %.2e\n", tgfb1_test$p.value))
  cat(sprintf("Mean TGFB1 - Hyperexpanded: %.4f, Other: %.4f\n",
              mean(low_hyperexp), mean(low_other)))
 }
 
 # Plot
 p_tgfb1 <- ggplot(tgfb1_data %>% filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize)),
                   aes(x = cloneSize, y = TGFB1, fill = BAFFreceptorStatus)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "TGFB1 Expression by Clonality",
       subtitle = "Do hyperexpanded LOW cells produce TGFB1?",
       x = "Clone Size", y = "TGFB1 Expression",
       fill = "BAFF Receptor") +
  scale_fill_manual(values = c("Low" = "#E41A1C", "High" = "#377EB8"))
 
 ggsave("tgfb1_by_clonality.pdf", p_tgfb1, width = 12, height = 6)
 cat("\nSaved: tgfb1_by_clonality.pdf\n")
}

# ============================================================================
# 3.2: MIF, BMP6, WNT Expression by Clonality
# ============================================================================

cat("\n3.2: Effector molecule expression by clonality\n")
cat(rep("-", 80), "\n\n")

# Key effector molecules
effector_genes <- c("MIF", "BMP6", "WNT10A", "WNT5B", "VEGFB", "IGF1")
available_effectors <- effector_genes[effector_genes %in% rownames(seuratObj)]

cat("Available effector genes:\n")
print(available_effectors)
cat("\n")

if (length(available_effectors) > 0) {
 
 # Extract data
 effector_data <- FetchData(seuratObj,
                            vars = c(available_effectors, "BAFFreceptorStatus", "cloneSize"))
 
 # Summary by clonality
 effector_summary <- effector_data %>%
  filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize)) %>%
  group_by(BAFFreceptorStatus, cloneSize) %>%
  summarise(across(all_of(available_effectors),
                   list(mean = mean, 
                        median = median,
                        pct = ~sum(. > 0) / n() * 100)),
            n = n(),
            .groups = "drop")
 
 write.csv(effector_summary, "effector_molecules_by_clonality.csv", row.names = FALSE)
 cat("Saved: effector_molecules_by_clonality.csv\n\n")
 
 # Key question: Are single HIGH cells the main producers?
 cat("Comparing effector production: Single vs Clonal HIGH receptor cells\n\n")
 
 high_single <- effector_data %>%
  filter(BAFFreceptorStatus == "High",
         cloneSize == "Single (0 < X <= 1)")
 
 high_clonal <- effector_data %>%
  filter(BAFFreceptorStatus == "High",
         cloneSize != "Single (0 < X <= 1)")
 
 for (gene in available_effectors) {
  cat(sprintf("\n%s:\n", gene))
  cat(sprintf("  Single HIGH cells: mean = %.4f, %% expressing = %.1f%%\n",
              mean(high_single[[gene]]),
              sum(high_single[[gene]] > 0) / nrow(high_single) * 100))
  
  if (nrow(high_clonal) > 0) {
   cat(sprintf("  Clonal HIGH cells: mean = %.4f, %% expressing = %.1f%%\n",
               mean(high_clonal[[gene]]),
               sum(high_clonal[[gene]] > 0) / nrow(high_clonal) * 100))
   
   # Test
   test <- wilcox.test(high_single[[gene]], high_clonal[[gene]])
   cat(sprintf("  p-value: %.2e\n", test$p.value))
  }
 }
 
 # Heatmap of effector expression by clonality
 effector_matrix <- effector_summary %>%
  select(BAFFreceptorStatus, cloneSize, ends_with("_mean")) %>%
  unite("group", BAFFreceptorStatus, cloneSize, sep = "_") %>%
  column_to_rownames("group") %>%
  as.matrix()
 
 # Clean column names
 colnames(effector_matrix) <- gsub("_mean$", "", colnames(effector_matrix))
 
 pdf("effector_molecules_clonality_heatmap.pdf", width = 10, height = 8)
 
 Heatmap(t(effector_matrix),
         name = "Mean\nExpression",
         col = colorRamp2(c(0, 0.1, 0.5, 1), c("white", "lightblue", "orange", "red")),
         cluster_rows = TRUE,
         cluster_columns = TRUE,
         row_names_side = "left",
         column_names_rot = 45,
         column_title = "Effector Molecule Expression by Clonality and Receptor Status",
         heatmap_legend_param = list(
          title = "Mean Expression"
         ))
 
 dev.off()
 cat("\nSaved: effector_molecules_clonality_heatmap.pdf\n")
}

# ============================================================================
# 3.3: Proportion of Producers by Clone Size
# ============================================================================

cat("\n3.3: What proportion of each clone size produces key effectors?\n")
cat(rep("-", 80), "\n\n")

if ("MIF" %in% rownames(seuratObj) & "BMP6" %in% rownames(seuratObj)) {
 
 producer_analysis <- FetchData(seuratObj,
                                vars = c("MIF", "BMP6", "TGFB1", 
                                         "BAFFreceptorStatus", "cloneSize")) %>%
  filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize)) %>%
  mutate(
   MIF_producer = MIF > 0,
   BMP6_producer = BMP6 > 0,
   TGFB1_producer = TGFB1 > 0
  ) %>%
  group_by(BAFFreceptorStatus, cloneSize) %>%
  summarise(
   n_cells = n(),
   pct_MIF = sum(MIF_producer) / n() * 100,
   pct_BMP6 = sum(BMP6_producer) / n() * 100,
   pct_TGFB1 = sum(TGFB1_producer) / n() * 100,
   .groups = "drop"
  )
 
 cat("Percentage of cells producing key molecules:\n")
 print(producer_analysis)
 
 # Plot
 producer_plot_data <- producer_analysis %>%
  pivot_longer(cols = starts_with("pct_"),
               names_to = "molecule",
               values_to = "percentage") %>%
  mutate(molecule = gsub("pct_", "", molecule))
 
 p_producers <- ggplot(producer_plot_data,
                       aes(x = cloneSize, y = percentage, fill = molecule)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~BAFFreceptorStatus, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Effector Molecule Production by Clonality",
       x = "Clone Size", y = "% Cells Expressing",
       fill = "Molecule") +
  scale_fill_brewer(palette = "Set2")
 
 ggsave("producer_frequency_by_clonality.pdf", p_producers, width = 12, height = 10)
 cat("\nSaved: producer_frequency_by_clonality.pdf\n")
}

# ============================================================================
# 3.4: Per-Clone Analysis (if clonotype IDs available)
# ============================================================================

cat("\n3.4: Per-clone effector molecule analysis\n")
cat(rep("-", 80), "\n\n")

if ("CTaa" %in% colnames(seuratObj@meta.data)) {
 
 # Get expression for all clones
 clone_expr <- FetchData(seuratObj,
                         vars = c("MIF", "BMP6", "WNT10A", "TGFB1",
                                  "BAFFreceptorStatus", "cloneSize", "CTaa")) %>%
  filter(!is.na(CTaa), CTaa != "", !is.na(BAFFreceptorStatus))
 
 # Calculate mean expression per clone
 per_clone_summary <- clone_expr %>%
  group_by(CTaa, BAFFreceptorStatus, cloneSize) %>%
  summarise(
   clone_size_n = n(),
   mean_MIF = mean(MIF),
   mean_BMP6 = mean(BMP6),
   mean_WNT10A = mean(WNT10A),
   mean_TGFB1 = mean(TGFB1),
   pct_MIF = sum(MIF > 0) / n() * 100,
   pct_BMP6 = sum(BMP6 > 0) / n() * 100,
   .groups = "drop"
  ) %>%
  arrange(desc(clone_size_n))
 
 cat("Top 20 clones by size:\n")
 print(per_clone_summary %>% head(20))
 
 write.csv(per_clone_summary, "per_clone_effector_analysis.csv", row.names = FALSE)
 cat("\nSaved: per_clone_effector_analysis.csv\n")
 
 # Identify "super-producer" clones
 super_producers <- per_clone_summary %>%
  filter(clone_size_n >= 10,  # At least 10 cells
         (pct_MIF > 50 | pct_BMP6 > 30))
 
 if (nrow(super_producers) > 0) {
  cat("\n'Super-producer' clones (>50% MIF+ or >30% BMP6+):\n")
  print(super_producers)
  
  write.csv(super_producers, "super_producer_clones.csv", row.names = FALSE)
  cat("\nSaved: super_producer_clones.csv\n")
 }
}

# ============================================================================
# SECTION 4: INTEGRATED SUMMARY
# ============================================================================

cat("\n\n")
cat("=" %R% 80, "\n")
cat("FINAL SUMMARY: KEY FINDINGS\n")
cat("=" %R% 80, "\n\n")

# Create comprehensive summary table
if (exists("effector_summary") & exists("tgfb1_summary")) {
 
 integrated_summary <- list(
  clone_distribution = table(seuratObj$BAFFreceptorStatus, seuratObj$cloneSize),
  
  key_findings = data.frame(
   Question = c(
    "Are LOW receptor cells more clonal?",
    "Do hyperexpanded LOW cells produce TGFB1?",
    "Are single HIGH cells main MIF producers?",
    "Are single HIGH cells main BMP6 producers?"
   ),
   Answer = c(
    sprintf("YES - %.1f%% in large/hyperexpanded clones vs %.1f%% for HIGH",
            (388+361)/1167*100, (22+105)/673*100),
    "Check tgfb1_summary table",
    "Check effector_summary table",
    "Check effector_summary table"
   )
  )
 )
 
 cat("\nKEY FINDINGS:\n")
 print(integrated_summary$key_findings)
 
 cat("\n\nCLONE SIZE DISTRIBUTION:\n")
 print(integrated_summary$clone_distribution)
 print(prop.table(integrated_summary$clone_distribution, margin = 1) * 100)
}

cat("\n\n")
cat("=" %R% 80, "\n")
cat("ANALYSIS COMPLETE\n")
cat("=" %R% 80, "\n\n")

cat("Generated files:\n")
cat("  - v_gene_usage_comparison.pdf\n")
cat("  - public_clones_shared.csv\n")
cat("  - aicda_expression_clonality.pdf\n")
cat("  - antigen_experience_markers_by_clonality.csv\n")
cat("  - tgfb1_by_clonality.pdf\n")
cat("  - effector_molecules_by_clonality.csv\n")
cat("  - effector_molecules_clonality_heatmap.pdf\n")
cat("  - producer_frequency_by_clonality.pdf\n")
cat("  - per_clone_effector_analysis.csv\n")
cat("  - super_producer_clones.csv (if applicable)\n")

cat("\n\nNEXT STEPS:\n")
cat("1. Review V gene usage patterns - are hyperexpanded clones using specific V genes?\n")
cat("2. Check AICDA expression - marker of active somatic hypermutation\n")
cat("3. Examine TGFB1 in hyperexpanded LOW cells - still regulatory?\n")
cat("4. Determine if single or clonal HIGH cells produce more MIF/BMP6\n")
cat("5. Identify any 'super-producer' clones\n")
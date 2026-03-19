
# Liver B Cell BCR Repertoire & Clonality Analysis
# This script:
# 1. Loads and combines BCR contig data across all liver samples
# 2. Adds BCR data to Seurat object
# 3. Analyses IGHV/IGLV gene usage by BAFF receptor status
# 4. Analyses BCR V(D)J gene usage, clonality, and somatic hypermutation
# 5. Correlates clonality with effector function

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(scRepertoire)
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
options(svglite.strict.text = FALSE)

# SECTION 1: LOAD AND COMBINE BCR CONTIG DATA

setwd('/data/home/hdx044/files/screpertoire/demux_contig/BCR')

# Healthy
s4  <- read.csv("GC-WL-10742-LIVER_LIVER_BCR_contig.csv")

# F0 (No fibrosis)
s8  <- read.csv("GC-WL-9961-LIVER_LIVER_BCR_contig.csv")
s12 <- read.csv("GC-WL-9999-LIVER_LIVER_BCR_contig.csv")
s16 <- read.csv("GC-WL-10113-2-LIVER_LIVER_BCR_contig.csv")

# F1 (Fibrosis)
s20 <- read.csv("GC-WL-9680-LIVER_LIVER_BCR_contig.csv")
s24 <- read.csv("GC-WL-10203-LIVER_LIVER_BCR_contig.csv")
s28 <- read.csv("GC-WL-10113-1-LIVER_LIVER_BCR_contig.csv")
s32 <- read.csv("GC-WL-10380-LIVER_LIVER_BCR_contig.csv")
s36 <- read.csv("GC-WL-10291-2-LIVER_LIVER_BCR_contig.csv")
s40 <- read.csv("GC-WL-10202-LIVER_LIVER_BCR_contig.csv")
s44 <- read.csv("GC-WL-10634-LIVER_LIVER_BCR_contig.csv")
s48 <- read.csv("GC-WL-10738-LIVER_LIVER_BCR_contig.csv")
s80 <- read.csv("GC-WL-11327-LIVER_LIVER_BCR_contig.csv")

# F2 (Fibrosis)
s52 <- read.csv("GC-WL-10205-LIVER_LIVER_BCR_contig.csv")
s56 <- read.csv("GC-WL-9991-LIVER_LIVER_BCR_contig.csv")
s60 <- read.csv("GC-WL-9932-LIVER_LIVER_BCR_contig.csv")
s64 <- read.csv("GC-WL-11040-LIVER_LIVER_BCR_contig.csv")

# F3 (Fibrosis)
s68 <- read.csv("GC-WL-11051-LIVER_LIVER_BCR_contig.csv")
s72 <- read.csv("GC-WL-11183-LIVER_LIVER_BCR_contig.csv")
s76 <- read.csv("GC-WL-11471-LIVER_LIVER_BCR_contig.csv")

liver_contig_list <- list(s4, s8, s12, s16, s20, s24, s28, s32, s36, s40,
                          s44, s48, s52, s56, s60, s64, s68, s72, s76, s80)

sample_names <- c("GC-WL-10742-LIVER", "GC-WL-9961-LIVER", "GC-WL-9999-LIVER",
                  "GC-WL-10113-2-LIVER", "GC-WL-9680-LIVER", "GC-WL-10203-LIVER",
                  "GC-WL-10113-1-LIVER", "GC-WL-10380-LIVER", "GC-WL-10291-2-LIVER",
                  "GC-WL-10202-LIVER", "GC-WL-10634-LIVER", "GC-WL-10738-LIVER",
                  "GC-WL-10205-LIVER", "GC-WL-9991-LIVER", "GC-WL-9932-LIVER",
                  "GC-WL-11040-LIVER", "GC-WL-11051-LIVER", "GC-WL-11183-LIVER",
                  "GC-WL-11471-LIVER", "GC-WL-11327-LIVER")

combined.BCR.liver <- combineBCR(liver_contig_list,
                                 samples = sample_names,
                                 removeNA = FALSE,
                                 removeMulti = FALSE,
                                 filterMulti = FALSE)

# Add disease type variable
combined.BCR.liver <- addVariable(combined.BCR.liver,
                                  variable.name = "Type",
                                  variables = c("Healthy", "No_fibrosis", "No_fibrosis",
                                                "No_fibrosis", "Fibrosis", "Fibrosis",
                                                "Fibrosis", "Fibrosis", "Fibrosis",
                                                "Fibrosis", "Fibrosis", "Fibrosis",
                                                "Fibrosis", "Fibrosis", "Fibrosis",
                                                "Fibrosis", "Fibrosis", "Fibrosis",
                                                "Fibrosis", "Fibrosis"))

# Add fibrosis stage variable
combined.BCR.liver <- addVariable(combined.BCR.liver,
                                  variable.name = "Fibrosis_Stage",
                                  variables = c("Healthy", "F0", "F0", "F0",
                                                "F1", "F1", "F1", "F1", "F1",
                                                "F1", "F1", "F1", "F2", "F2",
                                                "F2", "F2", "F3", "F3", "F3", "F1"))

setwd("/data/home/hdx044/files/screpertoire/BCell")

exportClones(combined.BCR.liver, write.file = TRUE, file.name = "allclonesliverBcell.csv")

# SECTION 2: ADD BCR DATA TO SEURAT OBJECT

seuratObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/liverBaffreceptorObj.rds')

# Fix barcode format
seurat_barcode_clean <- sub("_\\d+$", "", colnames(seuratObj))
seuratObj$barcode_matched <- paste0(seuratObj$orig.ident, "_", seurat_barcode_clean)
new_names <- seuratObj$barcode_matched
names(new_names) <- colnames(seuratObj)
seuratObj <- RenameCells(seuratObj, new.names = new_names)

# Combine BCR expression into Seurat object
seuratObj <- combineExpression(combined.BCR.liver,
                               seuratObj,
                               cloneCall = "gene",
                               group.by = "sample",
                               proportion = FALSE,
                               cloneSize = c(Single = 1, Small = 5, Medium = 20,
                                             Large = 100, Hyperexpanded = 500))

table(seuratObj$cloneSize, useNA = "ifany")

saveRDS(seuratObj, '/data/Blizard-AlazawiLab/rk/seurat/liverBaffreceptorObjBCR.rds')

# SECTION 3: JOIN BAFF STATUS INTO BCR DATA FRAME

combined.BCR.liver <- bind_rows(combined.BCR.liver)

combined.BCR.liver <- combined.BCR.liver %>%
  mutate(
    sample_short      = sub("_.*", "", barcode),
    cellbarcode_fixed = sub(".*_", "", barcode)
  )

seurat_clean          <- sub("_[0-9]+$", "", colnames(seuratObj))
seuratObj$cellbarcode <- sub(".*_", "", seurat_clean)

meta_lookup <- seuratObj@meta.data %>%
  mutate(cellbarcode = seuratObj$cellbarcode) %>%
  select(orig.ident, cellbarcode, BAFFreceptorStatus) %>%
  distinct()

combined.BCR.liver <- combined.BCR.liver %>%
  left_join(meta_lookup,
            by = c("sample_short" = "orig.ident",
                   "cellbarcode_fixed" = "cellbarcode"))

print(table(!is.na(combined.BCR.liver$BAFFreceptorStatus)))

# Filter to cells with BAFF status
combined.BCR.liver.baff <- combined.BCR.liver %>%
  filter(!is.na(BAFFreceptorStatus))

table(combined.BCR.liver.baff$BAFFreceptorStatus)

# Re-split for scRepertoire functions
bcr_list_baff <- split(combined.BCR.liver.baff, combined.BCR.liver.baff$sample_short)
bcr_list_baff <- bcr_list_baff[sapply(bcr_list_baff, nrow) > 0]

#### Clone quantification ####
colors <- c(
 "High" = "#E41A1C",     # red
 "Low"  = "#377EB8"      # blue
)
ggp <- clonalQuant(bcr_list_baff, cloneCall = "gene", group.by = "BAFFreceptorStatus", scale = TRUE) +
 scale_fill_manual(values = colors) +
 scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +  # Tick every 10
 ggtitle("B cell clones in liver") + 
 theme_classic(base_size = 12) +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.line = element_line(color = "black"),
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.key.size = unit(1.5, "lines")
 )

ggp

ggp@data

ggsave(
 filename = "uniqueBcellClonePercentLiver.png",
 plot = ggp,
 device = "png",
 path = "/data/home/hdx044/plots/BAFF",
 width = 4,           # in inches
 height = 5,          # in inches
 units = "in",
 dpi = 300            # 300–600 for print-quality
)

# SECTION 4: IGHV AND LIGHT CHAIN V-GENE USAGE

# --- 4.1: Heavy chain V-gene ---
bcr <- combined.BCR.liver %>%
  filter(!is.na(BAFFreceptorStatus)) %>%
  mutate(IGHV = if ("IGH" %in% colnames(.)) sub("\\..*$", "", IGH) else NA_character_) %>%
  filter(!is.na(IGHV), IGHV != "")

cnt <- bcr %>%
  count(BAFFreceptorStatus, IGHV, name = "n") %>%
  group_by(IGHV) %>% mutate(total = sum(n)) %>% ungroup() %>%
  mutate(IGHV = fct_reorder(IGHV, total, .desc = TRUE))

prop <- cnt %>%
  group_by(BAFFreceptorStatus) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

#### IGH-V GENE PLOT ####

ggp <- ggplot(prop, aes(x = IGHV, y = prop, fill = BAFFreceptorStatus)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "IGH V-gene usage by BAFF receptor status (proportion)",
       x = "IGHV", y = "Proportion of BCRs") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggp

colors <- c(
 "High" = "#E41A1C",     # red
 "Low"  = "#377EB8"      # blue
)

ggp <- ggplot(prop, aes(x = IGHV, y = prop, fill = BAFFreceptorStatus)) +
 geom_col(position = "dodge") +
 scale_fill_manual(values = colors) +
 scale_y_continuous(labels = scales::percent_format()) +
 labs(
  title = "IGH V-gene usage by BAFF receptor status (proportion)",
  x = "IGHV",
  y = "Proportion of BCRs"
 ) +
 theme_classic(base_size = 12) +
 theme(
  text = element_text(lineheight = 0.9),    # avoids overlap in PNG/SVG
  axis.text.x = element_text(angle = 60, hjust = 1),
  plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
  legend.position = "top",
  legend.title = element_blank(),
  plot.margin = margin(10, 10, 10, 10)
 )

ggp
# Note: IGHV1-69D, IGHV4-30-4, IGHV4-59 more used in LOW

ggsave(
 filename = "IGHVLiver.png",
 plot = ggp,
 device = "png",
 path = "/data/home/hdx044/plots/BAFF",
 width = 10,
 height = 6,
 units = "in",
 dpi = 300
)

#### IGKV IGLV GENE PLOT ####

bcr_light <- combined.BCR.liver %>%
  filter(!is.na(BAFFreceptorStatus)) %>%
  mutate(LightV = if ("IGLC" %in% colnames(.)) sub("\\..*$", "", IGLC) else NA_character_) %>%
  filter(!is.na(LightV), LightV != "")

cnt_light <- bcr_light %>%
  count(BAFFreceptorStatus, LightV, name = "n") %>%
  group_by(LightV) %>% mutate(total = sum(n)) %>% ungroup() %>%
  mutate(LightV = fct_reorder(LightV, total, .desc = TRUE))

prop_light <- cnt_light %>%
  group_by(BAFFreceptorStatus) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggp_light <- ggplot(prop_light, aes(x = LightV, y = prop, fill = BAFFreceptorStatus)) +
 geom_col(position = "dodge") +
 scale_fill_manual(values = colors_baff) +
 scale_y_continuous(labels = scales::percent_format()) +
 labs(
  title = "IG light-chain V-gene usage by BAFF receptor status (proportion)",
  x = "IGKV/IGLV",
  y = "Proportion of BCRs"
 ) +
 theme_classic(base_size = 12) +
 theme(
  text = element_text(lineheight = 0.9),
  axis.text.x = element_text(angle = 60, hjust = 1),
  plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
  legend.position = "top",
  legend.title = element_blank(),
  plot.margin = margin(10, 10, 10, 10)
 )

ggp_light
# Note: IGKV3-15, IGLV2-23, IGKV3-11 more used in LOW

ggsave(
 filename = "IGKVLiver.png",
 plot = ggp_light,
 device = "png",
 path = "/data/home/hdx044/plots/BAFF",
 width = 10,
 height = 6,
 units = "in",
 dpi = 300
)

# OUTPUT PATHS

plot_dir <- "/data/home/hdx044/plots/BAFF"
file_dir <- "/data/home/hdx044/files/BAFF"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)

# SECTION 5: BCR CLONALITY ANALYSIS

cat("=" , "\n")
cat("SECTION 5: BCR V(D)J GENE USAGE IN EXPANDED CLONES\n\n")

bcr_columns <- grep("^IG[HKL]|^CTaa|^CTnt|clone",
                    colnames(seuratObj@meta.data),
                    value = TRUE, ignore.case = TRUE)
cat("Available BCR-related columns:\n")
print(bcr_columns)

bcr_data <- seuratObj@meta.data %>%
 filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize))

hyperexpanded_low <- bcr_data %>%
 filter(BAFFreceptorStatus == "Low",
        cloneSize == "Hyperexpanded (100 < X <= 500)")
cat(sprintf("Hyperexpanded LOW receptor cells: %d\n", nrow(hyperexpanded_low)))

if ("IGH_v_gene" %in% colnames(bcr_data)) {
 
 v_gene_comparison <- bcr_data %>%
  filter(cloneSize %in% c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)")) %>%
  group_by(BAFFreceptorStatus, IGH_v_gene) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(BAFFreceptorStatus) %>%
  mutate(freq = n / sum(n) * 100) %>%
  arrange(BAFFreceptorStatus, desc(freq))
 
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
       x = "V Gene", y = "Frequency (%)", fill = "BAFF Receptor") +
  scale_fill_manual(values = c("Low" = "#E41A1C", "High" = "#377EB8"))
 
 ggsave(file.path(plot_dir, "v_gene_usage_comparison.pdf"), p_vgene, width = 12, height = 6)
}

if ("CTaa" %in% colnames(bcr_data)) {
 
 patient_col <- if ("patient" %in% colnames(bcr_data)) "patient" else "orig.ident"
 
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
 
 public_clones <- clone_sharing %>% filter(n_patients > 1)
 cat(sprintf("Total unique clonotypes: %d\n", nrow(clone_sharing)))
 cat(sprintf("Public clones (shared across patients): %d\n", nrow(public_clones)))
 
 if (nrow(public_clones) > 0) {
  write.csv(public_clones, file.path(file_dir, "public_clones_shared.csv"), row.names = FALSE)
  
  public_clone_receptor <- bcr_data %>%
   filter(CTaa %in% public_clones$CTaa) %>%
   group_by(BAFFreceptorStatus) %>%
   summarise(n = n(), .groups = "drop") %>%
   mutate(freq = n / sum(n) * 100)
  cat("\nReceptor status of cells in public clones:\n")
  print(public_clone_receptor)
 }
}

#### Clonal expansion plot ####

prop_data <- seuratObj@meta.data %>%
 filter(!is.na(cloneSize)) %>%
 group_by(BAFFreceptorStatus, cloneSize) %>%
 summarise(n = n(), .groups = 'drop')

ggp <- ggplot(prop_data, aes(x = BAFFreceptorStatus, y = n, fill = cloneSize)) +
 geom_bar(stat = "identity", position = "stack") +
 labs(y = "Cell count", title = "Clonal expansion") +
 theme_classic(base_size = 12) +
 theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.line = element_line(color = "black"),
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.key.size = unit(1.5, "lines")
 )

ggp

ggp@data

ggsave(file.path(plot_dir, "ExpandedClonotypeLiver.png"),
       plot = ggp, width = 6, height = 6, units = "in", dpi = 300)


#### Effector molecule production by clone size and BAFF receptor status ####

# Final gene list
effector_genes <- c(
 # ECM remodelling / collagen
 "COL4A3",    # Collagen type IV - basement membrane fibrosis
 "COL4A4",    # Collagen type IV - basement membrane fibrosis
 "PLOD1",     # Collagen crosslinking - fibrosis
 "CTHRC1",    # Collagen triple helix - HSC activation marker
 
 # Growth factors / fibrosis drivers
 "BMP6",      # Iron metabolism + HSC activation in MASH
 "IGF1",      # Hepatocyte survival + fibrosis
 "MDK",       # Midkine - promotes HSC activation + fibrosis
 "WNT5B",     # Wnt signalling + fibrosis
 "WNT4",      # Wnt signalling + fibrosis
 
 
 # Plasma cell / secretory markers
 "SDC1",      # CD138 - plasma cell + fibrosis niche
 "MZB1",      # Plasma cell ER function
 "XBP1",      # UPR / plasma cell differentiation
 
 # Signalling / transcription factors
 "ITGA8",     # Integrin alpha-8 - HSC activation
 "MERTK",     # TAM receptor - fibrosis + immune regulation
 "ESR1",      # Estrogen receptor - liver fibrosis
 "TGFBR3L",   # TGF-beta receptor - fibrosis signalling
 
 # Metabolism / stress (links to Fig 2 pathways)
 "PRDX4",     # Peroxiredoxin - ER oxidative stress + fibrosis
 "TXNDC5",    # Thioredoxin - ER stress + collagen folding
 "PHGDH",     # Serine synthesis - metabolic reprogramming
 
 # Keep for comparison (tells different story)
 "TGFB1",     # Canonical fibrosis - higher in LOW hyperexpanded
)

available_effectors <- effector_genes[effector_genes %in% rownames(seuratObj)]
cat("Available effector genes:\n")
print(available_effectors)

# Extract data
effector_data <- FetchData(seuratObj,
                           vars = c(available_effectors,
                                    "BAFFreceptorStatus", "cloneSize")) %>%
 filter(!is.na(BAFFreceptorStatus), !is.na(cloneSize))

# Compute % cells expressing per group
effector_pct <- effector_data %>%
 group_by(BAFFreceptorStatus, cloneSize) %>%
 summarise(across(all_of(available_effectors),
                  \(x) sum(x > 0) / n() * 100),
           n_cells = n(),
           .groups = "drop")

# Sanity check
cat("Groups and cell counts:\n")
print(effector_pct[, c("BAFFreceptorStatus", "cloneSize", "n_cells")])

# Pivot to long format
effector_long <- effector_pct %>%
 pivot_longer(cols = all_of(available_effectors),
              names_to = "gene",
              values_to = "pct_expressing") %>%
 mutate(
  # Shorten clone size labels for x axis
  cloneSize_short = case_when(
   grepl("Single",       cloneSize) ~ "Single",
   grepl("Small",        cloneSize) ~ "Small",
   grepl("Medium",       cloneSize) ~ "Medium",
   grepl("Large",        cloneSize) ~ "Large",
   grepl("Hyperexpanded",cloneSize) ~ "Hyper"
  ),
  cloneSize_short = factor(cloneSize_short,
                           levels = c("Single", "Small", "Medium",
                                      "Large", "Hyper"))
 )

# Plot
p_effector <- ggplot(effector_long,
                     aes(x = cloneSize_short, y = pct_expressing,
                         fill = BAFFreceptorStatus)) +
 geom_bar(stat = "identity", position = "dodge", width = 0.7) +
 facet_wrap(~gene, scales = "free_y", ncol = 5) +
 scale_fill_manual(values = c("High" = "#E41A1C", "Low" = "#377EB8")) +
 scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
 labs(
  title = "Effector molecule production by clone size and BAFF receptor status",
  x = "Clone size",
  y = "% Cells expressing",
  fill = "BAFF Receptor"
 ) +
 theme_classic(base_size = 12) +
 theme(
  strip.background = element_rect(fill = "grey90", color = NA),
  strip.text       = element_text(face = "bold", size = 11),
  axis.text.x      = element_text(angle = 45, hjust = 1),
  legend.position  = "top",
  legend.title     = element_text(size = 10),
  plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
  panel.spacing    = unit(1, "lines")
 )

# Save
ggsave(file.path(plot_dir, "effectorByCloneSize_BAFFstatus.png"),
       plot = p_effector, width = 10, height = 10, units = "in", dpi = 300)



























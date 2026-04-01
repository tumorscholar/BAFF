
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

#### load liver Baff receptor Obj BCR ####

seuratObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/liverBaffreceptorObjBCR.rds')

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

contingency_table <- table(seuratObj$BAFFreceptorStatus,
                           seuratObj$cloneSize)


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

#### IGH-V gene plot ####

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

#### IGKV IGLV gene plot ####

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


#### FNA BCR ####

setwd('/data/home/hdx044/files/screpertoire/demux_contig/BCR')

#Read TCR contig files

s3 = read.csv("GC-WL-10738-LIVER_LIVER_BCR_contig.csv")          # 10738 LIVER baseline
s4 = read.csv("GC-WL-11570-LIVER_BCR_contig.csv")                # 10738 LIVER followup

s6 = read.csv("GC-WL-10291-1-LIVER_BCR_contig.csv")              # 10291-1 LIVER baseline
s7 = read.csv("GC-WL-11303-LIVER_BCR_contig.csv")                # 10291-1 LIVER followup

s10 = read.csv("GC-WL-11040-LIVER_LIVER_BCR_contig.csv")         # 11040 LIVER baseline
s11 = read.csv("GC-WL-11816-LIVER_BCR_contig.csv")               # 11040 LIVER followup

s14 = read.csv("GC-WL-11183-LIVER_LIVER_BCR_contig.csv")         # 11183 LIVER baseline
s15 = read.csv("GC-WL-11937-LIVER_BCR_contig.csv")               # 11183 LIVER followup


# List of TCR data 
contig_list = list(s3,s4,s6,s7,s10,s11,s14,s15)

#Correct sample names (aligned with your real metadata)
sample_names <- c(
 "Baseline-10738-LIVER",
 "Followup-10738-LIVER",
 
 "Baseline-10291-1-LIVER",
 "Followup-10291-1-LIVER",
 
 "Baseline-11040-LIVER", 
 "Followup-11040-LIVER",
 
 "Baseline-11183-LIVER", 
 "Followup-11183-LIVER"
)

# Combine
combined.BCR <- combineBCR(
 contig_list,
 samples = sample_names,
 removeNA = FALSE,
 removeMulti = FALSE,
 filterMulti = FALSE
)

head(combined.BCR[[1]])

#### B cell count Baseline vs Followup ####

# Count B cells (unique barcodes) in each sample
bcell_counts <- sapply(combined.BCR, function(x) length(unique(x$barcode)))

# Put in a clean table
bcell_counts_df <- data.frame(
 Sample = names(combined.BCR),
 Bcell_Count = bcell_counts
)

bcell_counts_df

# Count B cells with productive IGH or IGL/IGK
bcell_productive <- sapply(combined.BCR, function(x) {
 sub <- x[!is.na(x$CTstrict), ]   # CTstrict = clonotype assignment
 length(unique(sub$barcode))
})

data.frame(
 Sample = names(combined.BCR),
 Productive_Bcells = bcell_productive
)

final_counts <- merge(bcell_counts_df,
                      data.frame(Sample = names(bcell_productive),
                                 Productive_Bcells = bcell_productive),
                      by = "Sample")

final_counts

total_cells <- c(
 "Baseline-10291-1-LIVER" = 1690,
 "Baseline-10738-LIVER"   = 3863,
 "Baseline-11040-LIVER"   = 1503,
 "Baseline-11183-LIVER"   = 2114,
 
 "Followup-10291-1-LIVER" = 3011,
 "Followup-10738-LIVER"   = 2348,
 "Followup-11040-LIVER"   = 2639,
 "Followup-11183-LIVER"   = 4366
)

# Add per-sample total cell counts
final_counts$Total_Cells <- total_cells[final_counts$Sample]

# Normalised B cells within each sample
final_counts$Bcell_fraction_within_sample <- final_counts$Bcell_Count / final_counts$Total_Cells

# Optional: Percentage
final_counts$Bcell_percent_within_sample <- final_counts$Bcell_fraction_within_sample * 100

#### B cell FNA final_counts ###
final_counts
#                   Sample Bcell_Count Productive_Bcells Total_Cells Bcell_fraction_within_sample Bcell_percent_within_sample
#1 Baseline-10291-1-LIVER          20                12        1690                   0.01183432                    1.183432
#2   Baseline-10738-LIVER         160               150        3863                   0.04141859                    4.141859
#3   Baseline-11040-LIVER         293               278        1503                   0.19494345                   19.494345
#4   Baseline-11183-LIVER         174               149        2114                   0.08230842                    8.230842
#5 Followup-10291-1-LIVER          53                44        3011                   0.01760213                    1.760213
#6   Followup-10738-LIVER          81                63        2348                   0.03449744                    3.449744
#7   Followup-11040-LIVER         196               147        2639                   0.07427056                    7.427056
#8   Followup-11183-LIVER         265               230        4366                   0.06069629                    6.069629

#### Clonal expansion analysis ####
# BCR Baseline vs Followup Analysis 

# Total B cells
total_bcells <- nrow(bind_rows(combined.BCR))
cat("Total number of B cells:", total_bcells, "\n")
#Total number of B cells: 1242

# Combine all samples
all_bcells <- dplyr::bind_rows(combined.BCR, .id = "Sample")

# Unique clonotypes
n_unique_clonotypes <- length(unique(all_bcells$CTaa))
cat("Number of unique clonotypes:", n_unique_clonotypes, "\n")
# Number of unique clonotypes: 867 

# Expanded vs singleton
clone_counts <- all_bcells %>%
 group_by(CTaa) %>%
 summarise(n_cells = n(), .groups = "drop")

expanded_clones  <- clone_counts %>% filter(n_cells > 1)
singleton_clones <- clone_counts %>% filter(n_cells == 1)

n_expanded_cells  <- sum(expanded_clones$n_cells)
n_singleton_cells <- sum(singleton_clones$n_cells)

cat("Expanded B cells:", n_expanded_cells, "\n")
# Expanded B cells: 419 

cat("Singleton B cells:", n_singleton_cells, "\n")
# Singleton B cells: 823

cat("Overall expansion percentage:", 
    round((n_expanded_cells / total_bcells) * 100, 2), "%\n")
# Overall expansion percentage: 33.74 %

# Add Patient and Timepoint metadata
for (i in seq_along(combined.BCR)) {
 sample_name <- names(combined.BCR)[i]
 
 if (grepl("10738", sample_name)) {
  combined.BCR[[i]]$Patient <- "10738"
 } else if (grepl("10291", sample_name)) {
  combined.BCR[[i]]$Patient <- "10291-1"
 } else if (grepl("11040", sample_name)) {
  combined.BCR[[i]]$Patient <- "11040"
 } else if (grepl("11183", sample_name)) {
  combined.BCR[[i]]$Patient <- "11183"
 }
 
 combined.BCR[[i]]$Timepoint <- ifelse(grepl("Baseline", sample_name),
                                       "Baseline", "Followup")
}

# Diversity analysis
div_shannon <- clonalDiversity(combined.BCR, cloneCall = "gene",
                               x.axis = "sample",
                               metric = "shannon", exportTable = TRUE)

div_inv_simpson <- clonalDiversity(combined.BCR, cloneCall = "gene",
                                   x.axis = "sample",
                                   metric = "inv.simpson", exportTable = TRUE)

prepare_div_data <- function(div_data, metric_name) {
 div_data$Patient <- ifelse(grepl("10738", div_data$sample), "10738",
                            ifelse(grepl("10291", div_data$sample), "10291-1",
                                   ifelse(grepl("11040", div_data$sample), "11040", "11183")))
 div_data$Timepoint <- ifelse(grepl("Baseline", div_data$sample),
                              "Baseline", "Followup")
 div_data$Patient   <- factor(div_data$Patient,
                              levels = c("10738", "10291-1", "11040", "11183"))
 div_data$Timepoint <- factor(div_data$Timepoint,
                              levels = c("Baseline", "Followup"))
 div_data$Metric    <- metric_name
 return(div_data)
}

shannon_data     <- prepare_div_data(div_shannon,     "Shannon")
inv_simpson_data <- prepare_div_data(div_inv_simpson, "Inverse Simpson")

timepoint_colors <- c("Baseline" = "#440154FF", "Followup" = "#FDE725FF")

p_shannon <- ggplot(shannon_data,
                    aes(x = Patient, y = value,
                        fill = Timepoint, color = Timepoint)) +
 geom_boxplot(alpha = 0.3, outlier.shape = NA) +
 geom_point(position = position_dodge(width = 0.75), size = 3) +
 geom_line(aes(group = Patient),
           position = position_dodge(width = 0.75),
           linetype = "dashed", alpha = 0.4, color = "black") +
 scale_fill_manual(values  = timepoint_colors) +
 scale_color_manual(values = timepoint_colors) +
 theme_classic() +
 theme(axis.text.x    = element_text(size = 10, face = "bold"),
       legend.position = "none") +
 labs(title = "Shannon Index", x = "Patient", y = "Shannon")

p_inv_simpson <- ggplot(inv_simpson_data,
                        aes(x = Patient, y = value,
                            fill = Timepoint, color = Timepoint)) +
 geom_boxplot(alpha = 0.3, outlier.shape = NA) +
 geom_point(position = position_dodge(width = 0.75), size = 3) +
 geom_line(aes(group = Patient),
           position = position_dodge(width = 0.75),
           linetype = "dashed", alpha = 0.4, color = "black") +
 scale_fill_manual(values  = timepoint_colors) +
 scale_color_manual(values = timepoint_colors) +
 theme_classic() +
 theme(axis.text.x     = element_text(size = 10, face = "bold"),
       legend.position = "bottom") +
 labs(title = "Inverse Simpson Index",
      x = "Patient", y = "Inv. Simpson",
      fill = "Timepoint", color = "Timepoint")

combined_plot <- (p_shannon | p_inv_simpson) +
 plot_annotation(
  title = "BCR Diversity Metrics - LIVER Baseline vs Followup",
  theme = theme(plot.title = element_text(face = "bold",
                                          size = 16, hjust = 0.5))
 )

print(combined_plot)

ggsave("/data/home/hdx044/plots/screpertoire/liver/FNA/BCR/BCR_Diversity_liver_BaselineFollowup.png",
       combined_plot, width = 10, height = 5, dpi = 300)


# Top 100 clone dominance
top_clone_analysis <- data.frame(
 Sample     = character(),
 Timepoint  = character(),
 Patient    = character(),
 Top100_Percent = numeric(),
 stringsAsFactors = FALSE
)

for (i in seq_along(combined.BCR)) {
 sample_data <- combined.BCR[[i]]
 sample_name <- names(combined.BCR)[i]
 
 clone_freq        <- table(sample_data$CTstrict)
 clone_freq_sorted <- sort(clone_freq, decreasing = TRUE)
 
 total_cells <- sum(clone_freq)
 top100 <- (sum(clone_freq_sorted[1:min(100, length(clone_freq_sorted))]) /
             total_cells) * 100
 
 timepoint <- ifelse(grepl("Baseline", sample_name), "Baseline", "Followup")
 patient   <- gsub("(Baseline-|Followup-)(.+)(-LIVER)", "\\2", sample_name)
 
 top_clone_analysis <- rbind(top_clone_analysis, data.frame(
  Sample        = sample_name,
  Timepoint     = timepoint,
  Patient       = patient,
  Top100_Percent = round(top100, 2)
 ))
}

cat("\n=== TOP 100 CLONE FREQUENCY BY SAMPLE ===\n")
print(top_clone_analysis)

# Summary
top100_comparison <- top_clone_analysis %>%
 group_by(Timepoint) %>%
 summarise(
  Mean_Top100 = round(mean(Top100_Percent), 2),
  SD_Top100   = round(sd(Top100_Percent), 2),
  n_samples   = n(),
  .groups     = "drop"
 )
cat("\n=== TOP 100 COMPARISON: Baseline vs Followup ===\n")
print(top100_comparison)

# Change per donor
top100_change <- top_clone_analysis %>%
 select(Patient, Timepoint, Top100_Percent) %>%
 pivot_wider(names_from = Timepoint, values_from = Top100_Percent) %>%
 mutate(
  Change         = Followup - Baseline,
  Percent_Change = round((Followup - Baseline) / Baseline * 100, 1),
  Direction      = case_when(
   Change < -2 ~ "Improved ↓",
   Change > 2  ~ "Worsened ↑",
   TRUE        ~ "Stable"
  )
 ) %>%
 filter(!is.na(Baseline) & !is.na(Followup))

cat("\n=== TOP 100 CHANGE BY DONOR ===\n")
print(top100_change)

#### B cell properties Baseline vs Followup ####

seuratObj <- readRDS("/data/Blizard-AlazawiLab/rk/seurat/SeuObjFNA_LIVER_processed.rds")

# Fix barcode format
cell.barcodes <- rownames(seuratObj[[]])

# removing the _1 at the end of the barcodes (adjust regex if your suffix differs)
cell.barcodes <- stringr::str_split(cell.barcodes, "_", simplify = TRUE)[,1]

# adding the prefix of the orig.ident to the barcodes, assuming that is the sample IDs
cell.barcodes <- paste0(seuratObj$orig.ident, "_", cell.barcodes)
seuratObj <- RenameCells(seuratObj, new.names = cell.barcodes)

prefix_vector <- c(
 "Baseline-10738-LIVER"   = "GC-WL-10738-LIVER",
 "Followup-10738-LIVER"   = "GC-WL-11570-LIVER",
 
 "Baseline-10291-1-LIVER" = "GC-WL-10291-1-LIVER",
 "Followup-10291-1-LIVER" = "GC-WL-11303-LIVER",
 
 "Baseline-11040-LIVER"   = "GC-WL-11040-LIVER",
 "Followup-11040-LIVER"   = "GC-WL-11816-LIVER",
 
 "Baseline-11183-LIVER"   = "GC-WL-11183-LIVER",
 "Followup-11183-LIVER"   = "GC-WL-11937-LIVER"
)

for (i in seq_along(combined.BCR)) {
 
 prefix <- prefix_vector[names(combined.BCR)[i]]
 
 combined.BCR[[i]]$pure <- sub(".*_", "", combined.BCR[[i]]$barcode)
 
 combined.BCR[[i]]$barcode <- paste0(prefix, "_", combined.BCR[[i]]$pure)
}

# Combine BCR expression into Seurat object
seuratObj <- combineExpression(combined.BCR,
                               seuratObj,
                               cloneCall = "gene",
                               group.by = "sample",
                               proportion = FALSE,
                               cloneSize = c(Single = 1, Small = 5, Medium = 20,
                                             Large = 100, Hyperexpanded = 500))

table(seuratObj$cloneSize, useNA = "ifany")

# subset BCR cells

bcr_cells <- seuratObj[, !is.na(seuratObj$cloneSize)]
bcr_cells
ncol(bcr_cells)   # should be 1027

saveRDS(bcr_cells, '/data/Blizard-AlazawiLab/rk/seurat/FNAliverBCR.rds')

# Memory B‑cell gene vector
memoryB_genes <- c(
 # Core memory B-cell markers
 "CD27", "CR2", "CD24", "CD83",
 
 # Naive/unswitched follicular B-cell markers
 "MS4A1", "CD79A", "CD79B", "CD19",
 "IGHM", "IGHD", "CXCR5",
 
 # BCR signaling / survival
 "POU2AF1", "SPIB", "BANK1", "BACH2", "BLK",
 "FCRL1", "FCRLA", "FCRL2",
 
 # B-cell regulation
 "BCL11A", "HVCN1", "RALGPS2", "OSBPL10"
)

bcr_cells <- AddModuleScore(
 bcr_cells,
 features = list(inflamB_genes),
 name = "InflamB"
)

bcr_cells <- AddModuleScore(
 bcr_cells,
 features = list(memoryB_genes),
 name = "MemoryB"
)

bcr_cells$VisitLabel <- factor(
 bcr_cells$VisitLabel,
 levels = c(
  "baseline 10291-1",
  "followup 10291-1",
  "baseline 10738",
  "followup 10738",
  "baseline 11040",
  "followup 11040",
  "baseline 11183",
  "followup 11183"
 )
)

VlnPlot2(bcr_cells, features = "MemoryB1", group.by = "VisitLabel")






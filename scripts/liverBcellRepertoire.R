library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(scRepertoire)
library(dplyr)
library(SeuratExtend)
library(stringr)
library(purrr)

#Read files
setwd('/data/home/hdx044/files/screpertoire/demux_contig/BCR')

# Healthy sample
s4 = read.csv("GC-WL-10742-LIVER_LIVER_BCR_contig.csv")

# F0 samples (No fibrosis)
s8 = read.csv("GC-WL-9961-LIVER_LIVER_BCR_contig.csv")
s12 = read.csv("GC-WL-9999-LIVER_LIVER_BCR_contig.csv")
s16 = read.csv("GC-WL-10113-2-LIVER_LIVER_BCR_contig.csv")

# F1 samples (Fibrosis)
s20 = read.csv("GC-WL-9680-LIVER_LIVER_BCR_contig.csv")
s24 = read.csv("GC-WL-10203-LIVER_LIVER_BCR_contig.csv")
s28 = read.csv("GC-WL-10113-1-LIVER_LIVER_BCR_contig.csv")
s32 = read.csv("GC-WL-10380-LIVER_LIVER_BCR_contig.csv")
s36 = read.csv("GC-WL-10291-2-LIVER_LIVER_BCR_contig.csv")
s40 = read.csv("GC-WL-10202-LIVER_LIVER_BCR_contig.csv")
s44 = read.csv("GC-WL-10634-LIVER_LIVER_BCR_contig.csv")
s48 = read.csv("GC-WL-10738-LIVER_LIVER_BCR_contig.csv")
s80 = read.csv("GC-WL-11327-LIVER_LIVER_BCR_contig.csv")

# F2 samples (Fibrosis)
s52 = read.csv("GC-WL-10205-LIVER_LIVER_BCR_contig.csv")
s56 = read.csv("GC-WL-9991-LIVER_LIVER_BCR_contig.csv")
s60 = read.csv("GC-WL-9932-LIVER_LIVER_BCR_contig.csv")
s64 = read.csv("GC-WL-11040-LIVER_LIVER_BCR_contig.csv")

# F3 samples (Fibrosis)
s68 = read.csv("GC-WL-11051-LIVER_LIVER_BCR_contig.csv")
s72 = read.csv("GC-WL-11183-LIVER_LIVER_BCR_contig.csv")
s76 = read.csv("GC-WL-11471-LIVER_LIVER_BCR_contig.csv")

# Create list of liver samples only
liver_contig_list = list(s4, s8, s12, s16, s20, s24, s28, s32, s36, s40, 
                         s44, s48, s52, s56, s60, s64, s68, s72, s76, s80)

# Combine BCR data for liver samples
combined.BCR.liver = combineBCR(liver_contig_list, 
                                samples = c("GC-WL-10742-LIVER",
                                            "GC-WL-9961-LIVER",
                                            "GC-WL-9999-LIVER",
                                            "GC-WL-10113-2-LIVER",
                                            "GC-WL-9680-LIVER",
                                            "GC-WL-10203-LIVER",
                                            "GC-WL-10113-1-LIVER",
                                            "GC-WL-10380-LIVER",
                                            "GC-WL-10291-2-LIVER",
                                            "GC-WL-10202-LIVER",
                                            "GC-WL-10634-LIVER",
                                            "GC-WL-10738-LIVER",
                                            "GC-WL-10205-LIVER",
                                            "GC-WL-9991-LIVER",
                                            "GC-WL-9932-LIVER",
                                            "GC-WL-11040-LIVER",
                                            "GC-WL-11051-LIVER",
                                            "GC-WL-11183-LIVER",
                                            "GC-WL-11471-LIVER",
                                            "GC-WL-11327-LIVER"),
                                removeNA = FALSE, 
                                removeMulti = FALSE,
                                filterMulti = FALSE)

# Add disease type variable
combined.BCR.liver = addVariable(combined.BCR.liver, 
                                 variable.name = "Type", 
                                 variables = c("Healthy",           # 10742
                                               "No_fibrosis",       # 9961
                                               "No_fibrosis",       # 9999
                                               "No_fibrosis",       # 10113-2
                                               "Fibrosis",          # 9680
                                               "Fibrosis",          # 10203
                                               "Fibrosis",          # 10113-1
                                               "Fibrosis",          # 10380
                                               "Fibrosis",          # 10291-2
                                               "Fibrosis",          # 10202
                                               "Fibrosis",          # 10634
                                               "Fibrosis",          # 10738
                                               "Fibrosis",          # 10205
                                               "Fibrosis",          # 9991
                                               "Fibrosis",          # 9932
                                               "Fibrosis",          # 11040
                                               "Fibrosis",          # 11051
                                               "Fibrosis",          # 11183
                                               "Fibrosis",          # 11471
                                               "Fibrosis"))         # 11327

# Add fibrosis stage variable for more detailed classification
combined.BCR.liver = addVariable(combined.BCR.liver, 
                                 variable.name = "Fibrosis_Stage", 
                                 variables = c("Healthy",    # 10742
                                               "F0",          # 9961
                                               "F0",          # 9999
                                               "F0",          # 10113-2
                                               "F1",          # 9680
                                               "F1",          # 10203
                                               "F1",          # 10113-1
                                               "F1",          # 10380
                                               "F1",          # 10291-2
                                               "F1",          # 10202
                                               "F1",          # 10634
                                               "F1",          # 10738
                                               "F2",          # 10205
                                               "F2",          # 9991
                                               "F2",          # 9932
                                               "F2",          # 11040
                                               "F3",          # 11051
                                               "F3",          # 11183
                                               "F3",          # 11471
                                               "F1"))         # 11327


setwd ("/data/home/hdx044/files/screpertoire/BCell")

exportClones(combined.BCR.liver, 
             write.file = TRUE,
             file.name = "allclonesliverBcell.csv")

head(combined.BCR.liver)

clonalQuant(combined.BCR.liver, 
            cloneCall="gene", 
            group.by = "Type", 
            scale = TRUE)

clonalHomeostasis(combined.BCR.liver, 
                  group.by = "Type",
                  cloneCall = "gene")

#### Add BCR data to liver Obj ####
seuratObj <- readRDS('/data/Blizard-AlazawiLab/rk/seurat/liverBaffreceptorObj.rds')

# Fix barcode format: remove "_1" suffix and add sample prefix
seurat_barcode_clean <- sub("_\\d+$", "", colnames(seuratObj))
seuratObj$barcode_matched <- paste0(seuratObj$orig.ident, "_", seurat_barcode_clean)

# Rename Seurat cells to match BCR format
new_names <- seuratObj$barcode_matched
names(new_names) <- colnames(seuratObj)
seuratObj <- RenameCells(seuratObj, new.names = new_names)

# Combine BCR expression data
seuratObj <- combineExpression(combined.BCR.liver, 
                               seuratObj, 
                               cloneCall = "gene", 
                               group.by = "sample", 
                               proportion = FALSE, 
                               cloneSize = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

# Check results
table(seuratObj$cloneSize, useNA = "ifany")

Seurat::DimPlot(seuratObj, group.by = "cloneSize") +
 scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))

saveRDS(seuratObj, '/data/Blizard-AlazawiLab/rk/seurat/liverBaffreceptorObjBCR.rds')

# Plot
DimPlot(seuratObj, group.by = "cloneSize")

# Check correlation between BAFF receptor status and clonal expansion

# Cross-tabulation
table(seuratObj$BAFFreceptorStatus, seuratObj$cloneSize)

# Chi-square test
chisq.test(table(seuratObj$BAFFreceptorStatus[!is.na(seuratObj$cloneSize)], 
                 seuratObj$cloneSize[!is.na(seuratObj$cloneSize)]))

# plot
prop_data <- seuratObj@meta.data %>%
 filter(!is.na(cloneSize)) %>%
 group_by(BAFFreceptorStatus, cloneSize) %>%
 summarise(n = n(), .groups = 'drop')

ggplot(prop_data, aes(x = BAFFreceptorStatus, y = n, fill = cloneSize)) +
 geom_bar(stat = "identity", position = "stack") +
 labs(y = "Cell count", 
      title = "Clonal Expansion by BAFF Receptor Status") +
 theme_classic()

#### BCR REPERTOIRE ANALYSIS ####

combined.BCR.liver <- bind_rows(combined.BCR.liver)

# Build keys on BCR side
# barcode looks like: "GC-WL-11327-LIVER_TTAGGCAGTTGATTGC-1"
combined.BCR.liver <- combined.BCR.liver %>%
 mutate(
  sample_short      = sub("_.*", "", barcode),     # "GC-WL-11327-LIVER"
  cellbarcode_fixed = sub(".*_", "", barcode)      # "TTAGGCAGTTGATTGC-1"
 )

# Build keys on Seurat side (extract the *pure* barcode from colnames)
# Seurat colnames look like: "GC-WL-11327-LIVER_TTAGGCAGTTGATTGC-1_1"
seurat_clean              <- sub("_[0-9]+$", "", colnames(seuratObj))    # drop trailing _1 / _2
seuratObj$cellbarcode     <- sub(".*_", "", seurat_clean)                # keep only "TTAGGCAGTTGATTGC-1"

# Build lookup from Seurat metadata
meta_lookup <- seuratObj@meta.data %>%
 mutate(cellbarcode = seuratObj$cellbarcode) %>%
 select(orig.ident, cellbarcode, BAFFreceptorStatus) %>%
 distinct()

# Join BAFF into combined.BCR.liver by (sample, pure barcode)
combined.BCR.liver <- combined.BCR.liver %>%
 left_join(
  meta_lookup,
  by = c(
   "sample_short"      = "orig.ident",
   "cellbarcode_fixed" = "cellbarcode"
  )
 )

# Quick checks
print(table(!is.na(combined.BCR.liver$BAFFreceptorStatus)))
head(combined.BCR.liver[, c("sample_short", "cellbarcode_fixed", "barcode", "BAFFreceptorStatus")], 10)


combined.BCR.liver.baff <- combined.BCR.liver %>%
 filter(!is.na(BAFFreceptorStatus))

# quick check
table(combined.BCR.liver.baff$BAFFreceptorStatus)
nrow(combined.BCR.liver.baff)










library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(scRepertoire)
library(dplyr)
library(SeuratExtend)
library(stringr)
library(purrr)


setwd ('/data/home/hdx044/files/screpertoire/demux_contig/BCR')

#Read files
#Healthy sample
s1 = read.csv("GC-WL-10742-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s2 = read.csv("GC-WL-10742-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s3 = read.csv("GC-WL-10742-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s4 = read.csv("GC-WL-10742-LIVER_LIVER_BCR_contig.csv")
#F0 samples
s5 = read.csv("GC-WL-9961-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s6 = read.csv("GC-WL-9961-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s7 = read.csv("GC-WL-9961-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s8 = read.csv("GC-WL-9961-LIVER_LIVER_BCR_contig.csv")
s9 = read.csv("GC-WL-9999-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s10 = read.csv("GC-WL-9999-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s11 = read.csv("GC-WL-9999-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s12 = read.csv("GC-WL-9999-LIVER_LIVER_BCR_contig.csv")
s13 = read.csv("GC-WL-10113-2-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s14 = read.csv("GC-WL-10113-2-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s15 = read.csv("GC-WL-10113-2-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s16 = read.csv("GC-WL-10113-2-LIVER_LIVER_BCR_contig.csv")
#F1 samples
s17 = read.csv("GC-WL-9680-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s18 = read.csv("GC-WL-9680-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s19 = read.csv("GC-WL-9680-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s20 = read.csv("GC-WL-9680-LIVER_LIVER_BCR_contig.csv")
s21 = read.csv("GC-WL-10203-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s22 = read.csv("GC-WL-10203-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s23 = read.csv("GC-WL-10203-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s24 = read.csv("GC-WL-10203-LIVER_LIVER_BCR_contig.csv")
#s25 = read.csv("GC-WL-10113-1-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s26 = read.csv("GC-WL-10113-1-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s27 = read.csv("GC-WL-10113-1-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s28 = read.csv("GC-WL-10113-1-LIVER_LIVER_BCR_contig.csv")
s29 = read.csv("GC-WL-10380-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s30 = read.csv("GC-WL-10380-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s31 = read.csv("GC-WL-10380-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s32 = read.csv("GC-WL-10380-LIVER_LIVER_BCR_contig.csv")
#s33 = read.csv("GC-WL-10291-2-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s34 = read.csv("GC-WL-10291-2-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s35 = read.csv("GC-WL-10291-2-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s36 = read.csv("GC-WL-10291-2-LIVER_LIVER_BCR_contig.csv")
s37 = read.csv("GC-WL-10202-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s38 = read.csv("GC-WL-10202-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s39 = read.csv("GC-WL-10202-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s40 = read.csv("GC-WL-10202-LIVER_LIVER_BCR_contig.csv")
s41 = read.csv("GC-WL-10634-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s42 = read.csv("GC-WL-10634-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s43 = read.csv("GC-WL-10634-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s44 = read.csv("GC-WL-10634-LIVER_LIVER_BCR_contig.csv")
s45 = read.csv("GC-WL-10738-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s46 = read.csv("GC-WL-10738-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s47 = read.csv("GC-WL-10738-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s48 = read.csv("GC-WL-10738-LIVER_LIVER_BCR_contig.csv")
#F2 samples
s49 = read.csv("GC-WL-10205-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s50 = read.csv("GC-WL-10205-SAT-VAT-PBMC_VAT_BCR_contig.csv")
#s51 = read.csv("GC-WL-10205-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s52 = read.csv("GC-WL-10205-LIVER_LIVER_BCR_contig.csv")
s53 = read.csv("GC-WL-9991-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s54 = read.csv("GC-WL-9991-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s55 = read.csv("GC-WL-9991-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s56 = read.csv("GC-WL-9991-LIVER_LIVER_BCR_contig.csv")
s57 = read.csv("GC-WL-9932-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s58 = read.csv("GC-WL-9932-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s59 = read.csv("GC-WL-9932-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s60 = read.csv("GC-WL-9932-LIVER_LIVER_BCR_contig.csv")
s61 = read.csv("GC-WL-11040-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s62 = read.csv("GC-WL-11040-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s63 = read.csv("GC-WL-11040-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s64 = read.csv("GC-WL-11040-LIVER_LIVER_BCR_contig.csv")
#F3 samples
s65 = read.csv("GC-WL-11051-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s66 = read.csv("GC-WL-11051-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s67 = read.csv("GC-WL-11051-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s68 = read.csv("GC-WL-11051-LIVER_LIVER_BCR_contig.csv")
s69 = read.csv("GC-WL-11183-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s70 = read.csv("GC-WL-11183-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s71 = read.csv("GC-WL-11183-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s72 = read.csv("GC-WL-11183-LIVER_LIVER_BCR_contig.csv")
s73 = read.csv("GC-WL-11471-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s74 = read.csv("GC-WL-11471-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s75 = read.csv("GC-WL-11471-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s76 = read.csv("GC-WL-11471-LIVER_LIVER_BCR_contig.csv")

#F1 samples
s77 = read.csv("GC-WL-11327-SAT-VAT-PBMC_SAT_BCR_contig.csv")
s78 = read.csv("GC-WL-11327-SAT-VAT-PBMC_VAT_BCR_contig.csv")
s79 = read.csv("GC-WL-11327-SAT-VAT-PBMC_PBMC_BCR_contig.csv")
s80 = read.csv("GC-WL-11327-LIVER_LIVER_BCR_contig.csv")

#list
contig_list = list(s1,s2,s3,s4, s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s26,s27,s28,s29,s30,s31,s32,s34,s35,s36,s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47,s48,s49,s50,s52,s53,s54,s55,s56,s57,s58,s59,s60,s61,s62,s63,s64,s65,s66,s67,s68,s69,s70,s71,s72,s73,s74,s75,s76,s77,s78,s79,s80)

combined.BCR = combineBCR(contig_list, samples = c("GC-WL-10742-SAT-VAT-PBMC-SAT", "GC-WL-10742-SAT-VAT-PBMC-VAT", "GC-WL-10742-SAT-VAT-PBMC-PBMC", "GC-WL-10742-LIVER-LIVER",
                                                   "GC-WL-9961-SAT-VAT-PBMC-SAT", "GC-WL-9961-SAT-VAT-PBMC-VAT", "GC-WL-9961-SAT-VAT-PBMC-PBMC", "GC-WL-9961-LIVER-LIVER",
                                                   "GC-WL-9999-SAT-VAT-PBMC-SAT", "GC-WL-9999-SAT-VAT-PBMC-VAT", "GC-WL-9999-SAT-VAT-PBMC-PBMC", "GC-WL-9999-LIVER-LIVER",
                                                   "GC-WL-10113-2-SAT-VAT-PBMC-SAT", "GC-WL-10113-2-SAT-VAT-PBMC-VAT", "GC-WL-10113-2-SAT-VAT-PBMC-PBMC", "GC-WL-10113-2-LIVER-LIVER",
                                                   "GC-WL-9680-SAT-VAT-PBMC-SAT", "GC-WL-9680-SAT-VAT-PBMC-VAT", "GC-WL-9680-SAT-VAT-PBMC-PBMC", "GC-WL-9680-LIVER-LIVER",
                                                   "GC-WL-10203-SAT-VAT-PBMC-SAT", "GC-WL-10203-SAT-VAT-PBMC-VAT", "GC-WL-10203-SAT-VAT-PBMC-PBMC", "GC-WL-10203-LIVER-LIVER",
                                                   "GC-WL-10113-1-SAT-VAT-PBMC-VAT", "GC-WL-10113-1-SAT-VAT-PBMC-PBMC", "GC-WL-10113-1-LIVER-LIVER",
                                                   "GC-WL-10380-SAT-VAT-PBMC-SAT", "GC-WL-10380-SAT-VAT-PBMC-VAT", "GC-WL-10380-SAT-VAT-PBMC-PBMC", "GC-WL-10380-LIVER-LIVER",
                                                   "GC-WL-10291-2-SAT-VAT-PBMC-VAT", "GC-WL-10291-2-SAT-VAT-PBMC-PBMC", "GC-WL-10291-2-LIVER-LIVER",
                                                   "GC-WL-10202-SAT-VAT-PBMC-SAT", "GC-WL-10202-SAT-VAT-PBMC-VAT", "GC-WL-10202-SAT-VAT-PBMC-PBMC", "GC-WL-10202-LIVER-LIVER",
                                                   "GC-WL-10634-SAT-VAT-PBMC-SAT","GC-WL-10634-SAT-VAT-PBMC-VAT","GC-WL-10634-SAT-VAT-PBMC-PBMC", "GC-WL-10634-LIVER-LIVER", 
                                                   "GC-WL-10738-SAT-VAT-PBMC-SAT","GC-WL-10738-SAT-VAT-PBMC-VAT","GC-WL-10738-SAT-VAT-PBMC-PBMC","GC-WL-10738-LIVER-LIVER", 
                                                   "GC-WL-10205-SAT-VAT-PBMC-SAT", "GC-WL-10205-SAT-VAT-PBMC-VAT", "GC-WL-10205-LIVER-LIVER",
                                                   "GC-WL-9991-SAT-VAT-PBMC-SAT", "GC-WL-9991-SAT-VAT-PBMC-VAT", "GC-WL-9991-SAT-VAT-PBMC-PBMC", "GC-WL-9991-LIVER-LIVER",
                                                   "GC-WL-9932-SAT-VAT-PBMC-SAT", "GC-WL-9932-SAT-VAT-PBMC-VAT", "GC-WL-9932-SAT-VAT-PBMC-PBMC", "GC-WL-9932-LIVER-LIVER",
                                                   "GC-WL-11040-SAT-VAT-PBMC-SAT", "GC-WL-11040-SAT-VAT-PBMC-VAT", "GC-WL-11040-SAT-VAT-PBMC-PBMC", "GC-WL-11040-LIVER-LIVER",
                                                   "GC-WL-11051-SAT-VAT-PBMC-SAT", "GC-WL-11051-SAT-VAT-PBMC-VAT", "GC-WL-11051-SAT-VAT-PBMC-PBMC", "GC-WL-11051-LIVER-LIVER",
                                                   "GC-WL-11183-SAT-VAT-PBMC-SAT", "GC-WL-11183-SAT-VAT-PBMC-VAT", "GC-WL-11183-SAT-VAT-PBMC-PBMC","GC-WL-11183-LIVER-LIVER",
                                                   "GC-WL-11471-SAT-VAT-PBMC-SAT", "GC-WL-11471-SAT-VAT-PBMC-VAT", "GC-WL-11471-SAT-VAT-PBMC-PBMC","GC-WL-11471-LIVER-LIVER",
                                                   "GC-WL-11327-SAT-VAT-PBMC-SAT", "GC-WL-11327-SAT-VAT-PBMC-VAT", "GC-WL-11327-SAT-VAT-PBMC-PBMC","GC-WL-11327-LIVER-LIVER"),
                          removeNA = FALSE, 
                          removeMulti = FALSE,
                          filterMulti = FALSE)

combined.BCR = addVariable(combined.BCR, 
                           variable.name = "Type", 
                           variables = c("Healthy", "Healthy", "Healthy", "Healthy",
                                         "No_fibrosis","No_fibrosis","No_fibrosis","No_fibrosis",
                                         "No_fibrosis","No_fibrosis","No_fibrosis","No_fibrosis",
                                         "No_fibrosis","No_fibrosis","No_fibrosis","No_fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis",
                                         "Fibrosis","Fibrosis","Fibrosis","Fibrosis"))


#Tissue information

combined.BCR <- lapply(combined.BCR, function(df) {
 df %>%
  mutate(
   Tissue = str_extract(sample, "(SAT|VAT|PBMC|LIVER)$")
  )
})

table(unlist(lapply(combined.BCR, function(x) x$Tissue)), useNA = "ifany")

head(combined.BCR[[1]])

setwd ("/data/home/hdx044/files/screpertoire/BCell")

exportClones(combined.BCR, 
             write.file = TRUE,
             file.name = "allclonesBcell.csv")

# Count total number of B cells across all samples
total_Bcells <- sum(sapply(combined.BCR, nrow))
cat("Total number of B cells:", total_Bcells, "\n")
# Total number of B cells: 7686 

# Combine all samples into one data frame
all_Bcells <- dplyr::bind_rows(combined.BCR, .id = "Sample")

# Count unique clonotypes based on CDR3 AA sequence (clone identifier)
n_unique_clonotypes <- length(unique(all_Bcells$CTaa))
cat("Number of unique clonotypes:", n_unique_clonotypes, "\n")
# Number of unique clonotypes: 4729

# Count cells per clonotype
clone_counts <- all_Bcells %>%
 group_by(CTaa) %>%
 summarise(n_cells = n(), .groups = "drop")

# Identify expanded clonotypes (>1 cell)
hyperexpanded_clones <- clone_counts %>% filter(n_cells > 20)
expanded_clones <- clone_counts %>% filter(n_cells > 1)
singleton_clones <- clone_counts %>% filter(n_cells == 1)

# Calculate summary statistics
n_hyperexpanded_clonotypes <- nrow(hyperexpanded_clones)
n_hyperexpanded_cells <- sum(hyperexpanded_clones$n_cells)

n_expanded_clonotypes <- nrow(expanded_clones)
n_expanded_cells <- sum(expanded_clones$n_cells)

n_singleton_clonotypes <- nrow(singleton_clones)
n_singleton_cells <- sum(singleton_clones$n_cells)

#### Print clone counts ####
cat("HyperExpanded clonotypes:", n_hyperexpanded_clonotypes, "\n")
# HyperExpanded clonotypes: 45
cat("HyperExpanded B cells:", n_hyperexpanded_cells, "\n")
# HyperExpanded B cells: 2180
cat("Expanded clonotypes:", n_expanded_clonotypes, "\n")
# Expanded clonotypes: 364 
cat("Expanded B cells:", n_expanded_cells, "\n")
# Expanded B cells: 3321 
cat("Singleton clonotypes:", n_singleton_clonotypes, "\n")
# Singleton clonotypes: 4365
cat("Singleton T cells:", n_singleton_cells, "\n")
# Singleton T cells: 4365
# Sanity check (expanded + singleton should equal total)
cat("Check total cells (expanded + singleton):", n_expanded_cells + n_singleton_cells, "\n")

# Extract the clonotype sequences (CTaa) of expanded clones
expanded_clone_list <- all_Bcells %>%
 filter(CTaa %in% expanded_clones$CTaa) %>%   # expanded clones only
 filter(!is.na(CTaa)) %>%
 group_by(CTaa) %>%
 summarise(
  n_cells = n(),                              # true cell count
  Samples = paste(sort(unique(sample)), collapse = ", "),
  .groups = "drop"
 ) %>%
 arrange(desc(n_cells))

# Write to CSV
write_csv(expanded_clone_list, "expanded_Bcells_Clonotypes_364.csv")

cat("✅ File 'exexpanded_Bcells_Clonotypes_364.csv written successfully.\n")

bcell_counts <- all_Bcells %>%
 filter(!is.na(Tissue)) %>%   # remove rare NA entries
 group_by(sample, Type, Tissue) %>%
 summarise(
  n_B_cells = n(),
  .groups = "drop"
 ) %>%
 arrange(sample, Type, Tissue)

head(bcell_counts)

# Export as CSV
write.csv(bcell_counts, "BcellCounSummary.csv", row.names = FALSE)

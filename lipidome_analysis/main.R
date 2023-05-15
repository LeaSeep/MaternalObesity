# Lipidome Analysis ----
## Overview ----
# This script includes the lipidome analysis Kupffer-Cells (KC) in wt and Hif1a ko
# setting.
# Including statistical analysis and summarised overview.

## Input ----
# (all in data Folder):
#   - normalized Lipidome measurements (Excel-Sheet), for preprocessing details
#     see manuscript, raw Data on Metabolomics Workbench
#   - Row annotation based on raw data (initial data from sample metadata separated)
#   - sample Table as provided


## Output----
# (output to script's folder):
#   - statistical analysis results for desired comparisons 
#     (rds object - displayed in table-format in provided html file)
#   - Heatmap Lipids

# SetUp ----
setwd("lipidome_analysis")
if(!("pics" %in% list.dirs())){
  dir.create("pics", showWarnings = FALSE)
}
# source relevant custom functions from utils folder
library(pheatmap)
library(SummarizedExperiment)
source("../utils/doSigLFCHeatmap.R")
source("../utils/doPCA.R")
output_result_list <- list()

# wild type KC ----
colorTheme <- c("#606060","#c12c38","#e0775f","#f3b694","#fce2d0")
names(colorTheme) <- c("CD_CD_HFD","HFD_CD_CD",
                       "HFD_CD_HFD", "HFD_HFD_CD","HFD_HFD_HFD")

RawData <- read.csv("../data/DataMatrix_wtLipidomics_2023_01_30.csv",row.names = 1)
RowAnno <- read.csv("../data/RowAnno_wtLipidomics_2023_01_30.csv",row.names = 1)
sampleAnno <- read.csv("../data/sampleTableLipid_wtLipidomics_2023_01_30.csv",row.names = 1)


sampleAnno$MaternalDiet <- "MaternalObese"
sampleAnno$MaternalDiet[grepl("^CD",sampleAnno$diet)] <- "MaternalLean"

sampleAnno$diet <- factor(
  sampleAnno$diet, 
  levels = c("CD_CD_CD","CD_CD_HFD","HFD_CD_CD",
             "HFD_CD_HFD","HFD_HFD_CD","HFD_HFD_HFD"), 
  ordered = F)

Lipid_SumExp <- SummarizedExperiment(
  assays  = RawData,
  rowData = RowAnno,
  colData = sampleAnno
  )



# remove anything constant
Lipid_SumExp <- Lipid_SumExp[which(apply(as.data.frame(assay(Lipid_SumExp)),1,sd)!=0),]

rowData(Lipid_SumExp)$CLASS_2 <- gsub("\\(.*$","",rownames(Lipid_SumExp))
rowData(Lipid_SumExp)[grepl("TAG",rowData(Lipid_SumExp)$CLASS_2),"CLASS_2"] <- rowData(Lipid_SumExp)[grepl("TAG",rowData(Lipid_SumExp)$CLASS_2),"CLASS"]


output_result_list[["Lipid_SumExp_wt"]] <- Lipid_SumExp
## Statistical Analysis & Heatmap ----
# Take out all TAGODD
# Those are most liekly not from mouse origin
Lipid_SumExp_woTAGODD = Lipid_SumExp[!grepl("TAGODD",rowData(Lipid_SumExp)$NAME),]

# To easen the comparison of wt vs ko data, the union of detected classes
# accross both measurements where identified and taken

subsetClasses=c("CE","Cer","DAG","DiHexCer,HexCer","LPC","LPC-O","LPE",
                "MAG","PA","PC","PC-O","PE","PE-O","PG","PS","SM","TACG_unsat","TACG_sat")

Lipid_SumExp_subset = Lipid_SumExp_woTAGODD[which(rowData(Lipid_SumExp_woTAGODD)$CLASS_2 %in% subsetClasses),]
old = rowData(Lipid_SumExp_woTAGODD)
new = rowData(Lipid_SumExp_subset)
print(paste0(
  "Metabolites dropped to Class Selection: ",
  round(((nrow(old)-nrow(new))/nrow(old))*100,2),"% (",
  nrow(old)-nrow(new),
  " mets, out of ",nrow(old),")"
))

results_wt <- doSigLFCHeatmap(
  as.data.frame(assay(Lipid_SumExp_subset)),
  as.data.frame(rowData(Lipid_SumExp_subset)),
  as.data.frame(colData(Lipid_SumExp_subset)),
  summarise_by = "CLASS_2", # column of RowAnno
  LFC_between = "diet", # column of sampleAnno
  FC_ctrl = "CD_CD_CD", # one level from colum selected in 'LFC_between',
  givenFilename = paste0("pics/","Heatmap_wt_",Sys.Date(),".png"),
  subset = c("TACG_unsat","TACG_sat","DAG", "CE", "MAG"), # put to NULL to get entire Heatmap
  colorTheme = colorTheme # needs to be named
)

output_result_list[["Lipidomics_wt"]] <- results_wt

## PCA ----
# Take the entire data set without class selection
Lipid_SumExp_log2 <- Lipid_SumExp
# remove anything 

colorTheme <- c("#c6c6c6","#606060","#c12c38","#e0775f","#f3b694","#fce2d0")
assay(Lipid_SumExp_log2) <- log2(assay(Lipid_SumExp_log2)+1)

WT_PCA <- doPCA(
  Lipid_SumExp_log2,
  colorTheme = colorTheme,
  shapeVar = "MaternalDiet", # one of colnames in colData(dds)
  colorVar = "diet" # one of colnames in colData(dds)
)
ggsave(filename = paste0("pics/","PCA_wt_",Sys.Date(),".png"), plot=WT_PCA)
ggsave(filename = paste0("pics/","PCA_wt_",Sys.Date(),".svg"), plot=WT_PCA)


# Hif1a KC ----
colorTheme = c("#1f78b4","#b2df8a","#33a02c",
               "#fdbf6f","#ff7f00","#fb9a99","#e31a1c")

names(colorTheme) <-c("CDCDCD_ko","CDCDHFD_wt","CDCDHFD_ko",
                      "HFDCDCD_wt","HFDCDCD_ko","HFDHFDHFD_wt","HFDHFDHFD_ko")

RawData <- read.csv("../data/DataMatrix_Hif1aLipidomics_2023_02_06.csv",row.names = 1)
RowAnno <- read.csv("../data/RowAnno_Hif1aLipidomics_2023_02_06.csv",row.names = 1)
sampleAnno <- read.csv("../data/sampleTableLipid_Hif1aLipidomics_2023_02_06.csv",row.names = 1)

sampleAnno$merged <- factor(
  sampleAnno$merged, 
  levels = c("CDCDCD_wt","CDCDCD_ko","CDCDHFD_wt","CDCDHFD_ko",
             "HFDCDCD_wt","HFDCDCD_ko","HFDHFDHFD_wt","HFDHFDHFD_ko"),
  ordered = F)

# replace all 'special'+','/','-' characters with '.'
rownames(sampleAnno) <- gsub("\\+|-|/","\\.",rownames(sampleAnno))

Lipid_SumExp_ko <- SummarizedExperiment(
  assays  = RawData,
  rowData = RowAnno,
  colData = sampleAnno
)


# remove anything constant
Lipid_SumExp_ko <- Lipid_SumExp_ko[which(apply(as.data.frame(assay(Lipid_SumExp_ko)),1,sd)!=0),]


rowData(Lipid_SumExp_ko)$CLASS_2 <- gsub("\\(.*$","",rownames(Lipid_SumExp_ko))
rowData(Lipid_SumExp_ko)[grepl("TAG",rowData(Lipid_SumExp_ko)$CLASS_2),"CLASS_2"] <- rowData(Lipid_SumExp_ko)[grepl("TAG",rowData(Lipid_SumExp_ko)$CLASS_2),"CLASS"]

output_result_list[["Lipid_SumExp_ko"]] <- Lipid_SumExp

## Statistical Analysis  & Heatmap ----
# Take out all TAGODD
# Those are most liekly not from mouse origin
Lipid_SumExp_woTAGODD = Lipid_SumExp_ko[!grepl("TAGODD",rowData(Lipid_SumExp_ko)$NAME),]

# To easen the comparison of wt vs ko data, the union of detected classes
# accross both measurements where identified and taken

Lipid_SumExp_subset = Lipid_SumExp_woTAGODD[which(rowData(Lipid_SumExp_woTAGODD)$CLASS_2 %in% subsetClasses),]
old = rowData(Lipid_SumExp_woTAGODD)
new = rowData(Lipid_SumExp_subset)
print(paste0(
  "Metabolites dropped to Class Intersection: ",
  round(((nrow(old)-nrow(new))/nrow(old))*100,2),"% (",
  nrow(old)-nrow(new),
  " mets, out of ",nrow(old),")"
))

results_ko <- doSigLFCHeatmap(
  as.data.frame(assay(Lipid_SumExp_subset)),
  as.data.frame(rowData(Lipid_SumExp_subset)),
  as.data.frame(colData(Lipid_SumExp_subset)),
  summarise_by = "CLASS_2", # column of RowAnno
  LFC_between = "merged", # column of sampleAnno
  FC_ctrl = "CDCDCD_wt", # one level from colum selected in 'LFC_between',
  givenFilename = paste0("pics/","Heatmap_ko_",Sys.Date(),".png"),
  subset = c("TACG_unsat","TACG_sat","DAG", "CE", "MAG"), # put to NULL to get entire Heatmap
  colorTheme = colorTheme # needs to be named
)


output_result_list[["Lipidomics_ko"]] <- results_ko

## PCA ----

Lipid_SumExp_log2 <- Lipid_SumExp_ko
# ensure correct colors
colData(Lipid_SumExp_log2)$merged <- factor(
  colData(Lipid_SumExp_log2)$merged, 
  levels = c("CDCDCD_wt","CDCDCD_ko","CDCDHFD_wt","CDCDHFD_ko",
             "HFDCDCD_wt","HFDCDCD_ko","HFDHFDHFD_wt","HFDHFDHFD_ko"),
  ordered = F)

colorTheme <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c", 
                "#fdbf6f","#ff7f00","#fb9a99","#e31a1c")
assay(Lipid_SumExp_log2) <- log2(assay(Lipid_SumExp_log2)+1)

KO_PCA <- doPCA(
  Lipid_SumExp_log2,
  colorTheme = colorTheme,
  shapeVar = "Maternal.diet", # one of colnames in colData(dds)
  colorVar = "merged" # one of colnames in colData(dds)
)
ggsave(filename = paste0("pics/","PCA_ko_",Sys.Date(),".png"), plot=KO_PCA)
ggsave(filename = paste0("pics/","PCA_ko_",Sys.Date(),".svg"), plot=KO_PCA)

# Save Results ----
saveRDS(output_result_list,"Lipidomic_stats_result.rds")

setwd("..")

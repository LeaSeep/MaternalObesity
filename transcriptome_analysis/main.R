# Transcriptomics Analysis ----
## Overview ----
# This script includes the transcritpome analysis of Hepatocytes (HC) as well as
# Kupffer-Cells (KC) - both Hif1a ko as well as 
# Including DE, ORA and Ligand-Receptor analysis. The Input for CoCena is generated here.
# The CoCena analysis is done in the respective Rmd, to be found in the same folder

## Input ----
# (all in data Folder):
#   - dds Objects created with script 'CreatingDDS.R', saved in data folder
#   - row annotation as provided, produced with `biomaRt::getBM` using ENSEMBL
#     Mus Musculus - downloaded on 05.07.2022 
#   - ligand-receptor interaction matrix (from CellTalk), parsed into rds object

## Output----
# (output to script's folder):
#   - DE-results object (KF & HC) - displayed in table-format in provided html file
#   - PCA plots
#   - Over-representation (ORA) plots for desired contrasts (KEGG, GO, HALLMARK)
#   - ORA-result objects
#   - Ligand-Receptor plot
#   - CoCena Heatmap
#   - CoCena Cluster ORA plots

# SetUp ----

setwd("transcriptome_analysis")

library("clusterProfiler")
library("msigdbr")
library("circlize")
library("DESeq2")
library("ggplot2")
library("sva")
# source relevant custom functions from utils folder
source("../utils/preprocessing_dds.R")
source("../utils/doPCA.R")
source("../utils/doORA.R")
source("../utils/doLigandReceptorPlot.R")
source("../utils/doComBat.R")

output_result_list <- list()

dds_KC <- readRDS("../data/DESeq_Obj_KC.rds")
dds_HC <- readRDS("../data/DESeq_Obj_HC.rds")

# Settings for both
colorTheme = c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
               "#fdbf6f","#ff7f00","#fb9a99","#e31a1c")

dds_KC$Merged <- factor(dds_KC$Merged,
                        levels=c("wt_cdcdcd_KC","ko_cdcdcd_KC",
                                 "wt_cdcdhfd_KC","ko_cdcdhfd_KC",
                                 "wt_hfdcdcd_KC","ko_hfdcdcd_KC",
                                 "wt_hfdhfdhfd_KC","ko_hfdhfdhfd_KC"))
# DE-Analysis ----

## KC ----
#remove transcript idenitifier
rownames(dds_KC) <- gsub("\\..*","",rownames(dds_KC))
rownames(dds_KC) <- gsub("\\..*","",rownames(dds_KC))

# Note there where to outliers removed namely IDs ("12553","12555") due to their
# outlier-behaviour (check out PCA if not removed), 12568 was removed as their
# was a miss-specification of the genotype and origin of the sample is not a 100%
# clear (most likely ko)

removeOutliers <- T

if(removeOutliers){
dds_KC <- dds_KC[,-c(
    which(colnames(dds_KC) %in% c("12553","12555","12569"))
    )]
}

design(dds_KC)=~Merged
dds_KC_preFilter <- dds_KC
dds_KC <- preprocessing(dds_KC,10,
                        protCodingOnly=T,
                        removeConstRows=T,
                        filterPerSample=F)
de_seq_result_KC <- DESeq(dds_KC) 

res_KC <- results(de_seq_result_KC,
                  contrast =c("Merged","ko_hfdcdcd_KC","wt_hfdcdcd_KC"),
                  alpha = 0.1)
summary(res_KC)
output_result_list[["dds_KC"]] <- res_KC

### PCA ----
colnames(colData(de_seq_result_KC))
vst_de_seq_result_KC <- vst(de_seq_result_KC,blind=T)
vst_de_seq_result_KC$Merged <- factor(vst_de_seq_result_KC$Merged,
                        levels=c("wt_cdcdcd_KC","ko_cdcdcd_KC",
                                 "wt_cdcdhfd_KC","ko_cdcdhfd_KC",
                                 "wt_hfdcdcd_KC","ko_hfdcdcd_KC",
                                 "wt_hfdhfdhfd_KC","ko_hfdhfdhfd_KC"))
KC_PCA <- doPCA(
  vst_de_seq_result_KC,
  colorTheme = colorTheme,
  shapeVar = "Maternal.diet", # one of colnames in colData(dds)
  colorVar = "Merged" # one of colnames in colData(dds)
  )
ggsave(filename = paste0("KC_PCA_",Sys.Date(),".png"), plot=KC_PCA)

## HC ----
#remove transcript idenitifier
rownames(dds_HC) <- gsub("\\..*","",rownames(dds_HC))
rownames(dds_HC) <- gsub("\\..*","",rownames(dds_HC))

dds_HC_preFilter <- dds_HC
dds_HC <- preprocessing(dds_HC,10)
design(dds_HC)=~Merged
de_seq_result_HC <- DESeq(dds_HC) 
res_HC <- results(de_seq_result_HC,
                  contrast = c("Merged","ko_hfdcdcd_HC","wt_hfdcdcd_HC"),
                  alpha = 0.1)
summary(res_HC)
output_result_list[["dds_HC"]] <- res_HC
### PCA ----
colnames(colData(de_seq_result_HC))
vst_de_seq_result_HC <- vst(de_seq_result_HC,blind=T)
vst_de_seq_result_HC$Merged <- factor(vst_de_seq_result_HC$Merged,
                                      levels=c("wt_cdcdcd_HC","ko_cdcdcd_HC",
                                               "wt_cdcdhfd_HC","ko_cdcdhfd_HC",
                                               "wt_hfdcdcd_HC","ko_hfdcdcd_HC",
                                               "wt_hfdhfdhfd_HC","ko_hfdhfdhfd_HC"))
HC_PCA <- doPCA(
  vst_de_seq_result_HC,
  colorTheme = colorTheme,
  shapeVar = "Maternal.diet", # one of colnames in colData(dds)
  colorVar = "Merged" # one of colnames in colData(dds)
)
ggsave(filename = paste0("HC_PCA_",Sys.Date(),".png"), plot=HC_PCA)

# ORA - Analysis ----
#Universe HC
universe_entrez <- bitr(rownames(res_HC),
                        fromType="ENSEMBL",
                        toType="ENTREZID",
                        OrgDb="org.Mm.eg.db")$ENTREZID

HC_ORA_UP <- doOra(
  rownames(res_HC[which(res_HC$log2FoldChange>0 & res_HC$padj <0.1),]), # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=4,
  universe_entrez,
  filename = "ORA_DE/HC_UP" # will be ORA_[filename][type].png
  )

HC_ORA_DOWN <- doOra(
  rownames(res_HC[which(res_HC$log2FoldChange<0 & res_HC$padj <0.1),]), # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=4,
  universe_entrez,
  filename = "ORA_DE/HC_DOWN"
)

KC_ORA_UP <- doOra(
  rownames(res_KC[which(res_KC$log2FoldChange>0 & res_KC$padj <0.1),]), # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=4,
  universe_entrez,
  filename = "ORA_DE/KC_UP"
)
KC_ORA_DOWN <- doOra(
  rownames(res_KC[which(res_KC$log2FoldChange<0 & res_KC$padj <0.1),]), # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=6,
  universe_entrez,
  filename = "ORA_DE/KC_DOWN"
)

output_result_list[["ORA_HC"]] <- list(UP = HC_ORA_UP,
                                       DOWN = HC_ORA_DOWN
                                       )  

output_result_list[["ORA_KC"]] <- list(UP = KC_ORA_UP,
                                       DOWN = KC_ORA_DOWN
                                       ) 
  
# Ligand-Receptor Analysis ----
lr_mouse_table=read.table("../data/mouse_lr_pair_fromCellTalk.txt",header = T)
DE_LigandReceptor <- doLigandReceptorPlot(
  allDEGenes_ligand = res_KC[which(res_KC$padj < 0.1 & abs(res_KC$log2FoldChange)>2),],
  dds_obj_Ligand = dds_KC,
  allPresent_receptor = dds_HC,
  colorVar = "log2FoldChange",
  adjMatrix_LigandReceptor = lr_mouse_table
  )

# CoCena ----
## HC Co-expression analysis ----
# This was done closely following this: https://github.com/MarieOestreich/hCoCena/tree/main/showcase
# main and satellite showcase Notebooks
# The Input is generated in here and saved to data folder

# Note we apply an additional filter step, to remove low count genes accross
# all samples, precisely everything is kept, that
# as in 25% of samples >10 counts

#keep <- rowSums(counts(dds_HC) >= 10) >= ceiling(0.25*ncol(dds_HC))
#dds_HC_coCena <- dds_HC[keep,]
CoCena_Input_HC <- DESeq(dds_HC)

saveRDS(CoCena_Input_HC,"../data/CoCena_Input_HC.rds")

## KC Co-expression analysis ----

#keep <- rowSums(counts(dds_KC) >= 10) >= ceiling(0.25*ncol(dds_KC))
#dds_KC_coCena <- dds_KC[keep,]
CoCena_Input_KC <- DESeq(dds_KC)

saveRDS(CoCena_Input_KC,"../data/CoCena_Input_KC.rds")

# Save everything ----
saveRDS(output_result_list, file = "Transcriptomics_results.rds")

setwd("..")



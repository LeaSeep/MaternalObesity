# Transcriptomics Analysis ----
## Overview ----
# This script includes the transcriptome analysis of Hepatocytes (HC) as well as
# Kupffer-Cells (KC) - both Hif1a ko as well as wild type. 
# Including DE, ORA and Ligand-Receptor analysis. The Input for CoCena is generated here.
# The CoCena analysis is done in the respective Rmd, to be found in the same folder.
# The script is structured in sections, outline visible on the right

## Input ----
# (all in data Folder):
#   - dds Objects created with script 'CreatingDDS.R', saved in data folder
#   - row annotation as provided, produced with `biomaRt::getBM` using ENSEMBL
#     Mus Musculus - downloaded on 05.07.2022 
#   - ligand-receptor interaction matrix (from CellTalk), parsed into rds object

## Output----
# (output to script's folder):
#   - DE-results object (KC [wt and ko] & HC) - displayed in table-format in provided html file
#   - PCA plots
#   - Over-representation (ORA) plots for desired contrasts (KEGG, GO, HALLMARK)
#   - ORA-result objects
#   - Ligand-Receptor plot
#.  - CoCena Input


# SetUp ----

if(gsub("^.+/","",getwd())!="transcriptome_analysis"){
  setwd("transcriptome_analysis")
}

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

dds_KC_WT <- readRDS("../data/DESeq_Obj_KC_WT.rds")
dds_KC_WT_P0 <- readRDS("../data/DESeq_Obj_KC_WT_P0.rds")
dds_KC_WT_E14_5 <- readRDS("../data/DESeq_Obj_KC_WT_E14_5.rds")

# Settings for both
colorTheme = c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
               "#fdbf6f","#ff7f00","#fb9a99","#e31a1c")


colorTheme_wt <- c("#c6c6c6","#606060","#c32b38","#fee2d1","#8a4094","#ebdeec")
dds_KC$Merged <- factor(dds_KC$Merged,
                        levels=c("wt_cdcdcd_KC","ko_cdcdcd_KC",
                                 "wt_cdcdhfd_KC","ko_cdcdhfd_KC",
                                 "wt_hfdcdcd_KC","ko_hfdcdcd_KC",
                                 "wt_hfdhfdhfd_KC","ko_hfdhfdhfd_KC"))

# DE-Analysis ----

## KC ----
### Hif1a ----
#remove transcript idenitifier
rownames(dds_KC) <- gsub("\\..*","",rownames(dds_KC))
rownames(dds_KC) <- gsub("\\..*","",rownames(dds_KC))

# Note there where to outliers removed namely IDs ("12553","12555") due to their
# outlier-behaviour (check out PCA if not removed), 12569 was removed as their
# was a miss-specification of the genotype and origin of the sample is not a 100%
# clear (most likely ko)

removeOutliers <- F

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
                        filterPerSample=T)
de_seq_result_KC <- DESeq(dds_KC) 


res_KC <- results(de_seq_result_KC,
                  contrast =c("Merged","ko_hfdcdcd_KC","wt_hfdcdcd_KC"),
                  alpha = 0.1)
res_KC <- results(de_seq_result_KC,
                  contrast =c("Merged","ko_hfdcdcd_KC","wt_hfdcdcd_KC"),
                  alpha = 0.1)

test<-lfcShrink(de_seq_result_KC,contrast =c("Merged","ko_hfdcdcd_KC","wt_hfdcdcd_KC"),type="ashr",lfcThreshold=1)

summary(res_KC)
output_result_list[["dds_KC"]] <- res_KC

#### PCA ----
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
ggsave(filename = paste0("KC_PCA_",Sys.Date(),".png"), plot = KC_PCA)
ggsave(filename = paste0("KC_PCA_",Sys.Date(),".svg"), plot = KC_PCA)

### WT ----
#remove transcript idenitifier
rownames(dds_KC_WT) <- gsub("\\..*","",rownames(dds_KC_WT))
rownames(dds_KC_WT) <- gsub("\\..*","",rownames(dds_KC_WT))
dds_KC_WT$Condition <- gsub("H","h",dds_KC_WT$Condition)

dds_KC_WT$Condition <- as.factor(dds_KC_WT$Condition)

removeOutliers <-T
dim(dds_KC_WT)
if(removeOutliers){
  dds_KC_WT <- dds_KC_WT[,-c(
    which(colnames(dds_KC_WT) %in% c(
      "5819","5820","5821","5822","5823","5824","5825","5953"
      ))
  )]
}
dim(dds_KC_WT)

design(dds_KC_WT)=~Condition
dds_KC_preFilter <- dds_KC_WT
dds_KC_WT <- preprocessing(dds_KC_WT,10,
                        protCodingOnly=T,
                        removeConstRows=T,
                        filterPerSample=T)
de_seq_result_KC_WT <- DESeq(dds_KC_WT) 

res_KC_wt <- results(de_seq_result_KC_WT,
                  contrast =c("Condition","hfdcdcd","cdcdcd"),
                  alpha = 0.1)
summary(res_KC_wt)
output_result_list[["dds_KC_WT"]] <- res_KC_wt

#### PCA ----
colnames(colData(de_seq_result_KC_WT))
vst_de_seq_result_KC_wt <- vst(de_seq_result_KC_WT,blind=T)
vst_de_seq_result_KC_wt$Condition <- factor(vst_de_seq_result_KC_wt$Condition,
                                      levels=c("cdcdcd","cdcdhfd",
                                               "hfdcdcd","hfdhfdcd",
                                               "hfdcdhfd","hfdhfdhfd"),ordered = T)


KC_PCA <- doPCA(
  vst_de_seq_result_KC_wt,
  colorTheme = colorTheme_wt,
  shapeVar = "Maternal_diet", # one of colnames in colData(dds)
  colorVar = "Condition" # one of colnames in colData(dds)
)
ggsave(filename = paste0("KC_WT_PCA_",Sys.Date(),".png"), plot=KC_PCA)
ggsave(filename = paste0("KC_WT_PCA_",Sys.Date(),".svg"), plot=KC_PCA)


### WT P0 ----
#remove transcript idenitifier
rownames(dds_KC_WT_P0) <- gsub("\\..*","",rownames(dds_KC_WT_P0))
rownames(dds_KC_WT_P0) <- gsub("\\..*","",rownames(dds_KC_WT_P0))

dds_KC_WT_P0$Condition <- as.factor(dds_KC_WT_P0$Condition)
dds_KC_WT_P0$date <- as.factor(dds_KC_WT_P0$date)
table(dds_KC_WT_P0$date ,dds_KC_WT_P0$Condition)

dds_KC_WT_P0$date_batch <- "Day1"
dds_KC_WT_P0$date_batch[(which(dds_KC_WT_P0$Condition=="CD" & dds_KC_WT_P0$date=="20180820"))] <- "Day2" 
dds_KC_WT_P0$date_batch[(which(dds_KC_WT_P0$Condition=="HFD" & dds_KC_WT_P0$date=="20180814"))] <- "Day2" 
table(dds_KC_WT_P0$date_batch,dds_KC_WT_P0$Condition)


dds_KC_WT_P0$date_batch <- as.factor(dds_KC_WT_P0$date_batch)

dim(dds_KC_WT_P0)

design(dds_KC_WT_P0)=~Condition + date_batch
dds_KC_P0_preFilter <- dds_KC_WT_P0
dds_KC_WT_P0 <- preprocessing(dds_KC_WT_P0,10,
                           protCodingOnly=T,
                           removeConstRows=T,
                           filterPerSample=T)
de_seq_result_KC_WT_P0 <- DESeq(dds_KC_WT_P0) 
deg <- de_seq_result_KC_WT_P0 %>%
  select(SYMBOL, log2FoldChange, stat, pvalue,padj) %>% 
  filter(!is.na(stat)) %>% 
  column_to_rownames(var = "SYMBOL") %>%
  as.matrix()
deg <- as.data.frame(deg)
deg <- deg[!is.na(deg$log2FoldChange),]

deg <- deg[deg[,"pvalue"]<0.05,]
res_KC_wt_P0 <- results(de_seq_result_KC_WT_P0,
                     contrast =c("Condition","HFD","CD"),
                     alpha = 0.1)
summary(res_KC_wt_P0)
output_result_list[["dds_KC_WT_P0"]] <- res_KC_wt_P0

# save objects to use for TF - analysis
saveRDS(dds_KC_WT_P0,file = "../data/dds_KC_WT_P0.rds")
saveRDS(res_KC_wt_P0,file = "../data/res_KC_WT_P0.rds")

#### PCA ----
colnames(colData(de_seq_result_KC_WT_P0))

vst_de_seq_result_KC_wt_P0 <- vst(dds_KC_WT_P0,blind=T)

correctedObj <- doBatchCorrection(
  SumExp_obj = vst_de_seq_result_KC_wt_P0,
  design_factor = "Condition",
  batch_factor = "date_batch"
)

KC_PCA <- doPCA(
  correctedObj,
  colorTheme = colorTheme_wt,
  shapeVar = "Condition", # one of colnames in colData(dds)
  colorVar = "date", # one of colnames in colData(dds)
  xPC = "PC1",
  yPC = "PC2"
)
ggsave(filename = paste0("KC_WT_P0_PCA_PC1_2",Sys.Date(),".png"), plot=KC_PCA)
ggsave(filename = paste0("KC_WT_P0_PCA_PC1_2",Sys.Date(),".svg"), plot=KC_PCA)

### WT E14.5 ----
#remove transcript idenitifier
rownames(dds_KC_WT_E14_5) <- gsub("\\..*","",rownames(dds_KC_WT_E14_5))
rownames(dds_KC_WT_E14_5) <- gsub("\\..*","",rownames(dds_KC_WT_E14_5))

dds_KC_WT_E14_5$Condition <- as.factor(dds_KC_WT_E14_5$Condition)


dim(dds_KC_WT_E14_5)
dds_KC_WT_E14_5$date <- as.factor(dds_KC_WT_E14_5$date)
design(dds_KC_WT_E14_5)=~Condition + date
dds_KC_P0_preFilter <- dds_KC_WT_E14_5
dds_KC_WT_E14_5 <- preprocessing(dds_KC_WT_E14_5,10,
                              protCodingOnly=T,
                              removeConstRows=T,
                              filterPerSample=T)
de_seq_result_KC_WT_E14_5 <- DESeq(dds_KC_WT_E14_5) 

res_KC_wt_E14_5 <- results(de_seq_result_KC_WT_E14_5,
                        contrast =c("Condition","HFD","CD"),
                        alpha = 0.1)
summary(res_KC_wt_E14_5)
output_result_list[["dds_KC_WT_E14_5"]] <- res_KC_wt_E14_5

#### PCA ----
colnames(colData(dds_KC_WT_E14_5))
vst_de_seq_result_KC_wt_E14_5 <- vst(dds_KC_WT_E14_5,blind=T)
vst_de_seq_result_KC_wt_E14_5$Condition <- factor(vst_de_seq_result_KC_wt_E14_5$Condition,
                                               levels=c("CD","HFD"),ordered = T)


correctedObj <- doBatchCorrection(
  SumExp_obj = vst_de_seq_result_KC_wt_E14_5,
  design_factor = "Condition",
  batch_factor = "date"
)


KC_PCA <- doPCA(
  correctedObj,
  colorTheme = colorTheme_wt,
  shapeVar = "Condition", # one of colnames in colData(dds)
  colorVar = "date", # one of colnames in colData(dds)
  xPC = "PC1",
  yPC = "PC2"
)
ggsave(filename = paste0("KC_WT_E13_5_PCA_",Sys.Date(),".png"), plot=KC_PCA)
ggsave(filename = paste0("KC_WT_E13_5_PCA_",Sys.Date(),".svg"), plot=KC_PCA)



## HC ----
#remove transcript idenitifier
rownames(dds_HC) <- gsub("\\..*","",rownames(dds_HC))
rownames(dds_HC) <- gsub("\\..*","",rownames(dds_HC))

dds_HC_preFilter <- dds_HC
dds_HC <- preprocessing(dds_HC,
                        10,
                        protCodingOnly=T,
                        removeConstRows=T,
                        filterPerSample=T
                        )
design(dds_HC)=~Merged
de_seq_result_HC <- DESeq(dds_HC) 
res_HC <- results(de_seq_result_HC,
                  contrast = c("Merged","wt_hfdcdcd_HC","wt_cdcdcd_HC"), #c("Merged","ko_hfdcdcd_HC","wt_hfdcdcd_HC"),
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
ggsave(filename = paste0("HC_PCA_",Sys.Date(),".svg"), plot=HC_PCA)

### SingleGene Vis ----
selectedSymbols <- c("Lect2","Cxcl12","Ccl25","C3","Agt")

idx_selected <- which(rowData(vst_de_seq_result_HC)[,"SYMBOL"] %in% selectedSymbols)
data2plot <- as.data.frame(assay(vst_de_seq_result_HC[idx_selected,]))

#### Boxplots ----
data2plot$Gene <- rowData(vst_de_seq_result_HC)[rownames(data2plot),"SYMBOL"]
data2plot <- reshape2::melt(data2plot,value.name = "counts")
data2plot$group <- "none"
for(j in 1:length(data2plot$variable)){
  data2plot$group[j] <- as.character(colData(vst_de_seq_result_HC)[data2plot$variable[j],"Merged"])
}

data2plot$group <- factor(data2plot$group,
                          levels=c("wt_cdcdcd_HC","ko_cdcdcd_HC",
                                   "wt_cdcdhfd_HC","ko_cdcdhfd_HC",
                                   "wt_hfdcdcd_HC","ko_hfdcdcd_HC",
                                   "wt_hfdhfdhfd_HC","ko_hfdhfdhfd_HC"))

P_boxplots <- ggplot(data2plot, 
                     aes(y=counts,
                         x=group,
                         fill=group))+
  geom_boxplot()+
  geom_point(shape = 21,size=2,alpha=0.6)+
  scale_fill_manual(values = colorTheme,
                    name = "group")+
  xlab("Gene")+
  ylab("vst Counts")+
  theme_bw()+
  facet_wrap(Gene~.)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size=20))

ggsave(filename = paste0("HC_SelectedGenesBoxplot_",Sys.Date(),".png"), plot=P_boxplots)
ggsave(filename = paste0("HC_SelectedGenesBoxplot_",Sys.Date(),".svg"), plot=P_boxplots)

#### Heatmap ----
data2sum <- vst_de_seq_result_HC
print(length(selectedSymbols))
union_selected <-  union_selected <- rownames(rowData(CoCena_Input_HC)[rowData(CoCena_Input_HC)[,"SYMBOL"] %in% selectedSymbols,])

subsetData <- data2sum[rownames(data2sum) %in% union_selected,]
subsetData_df <- as.data.frame(t(assay(subsetData)))
subsetData_df <- cbind(subsetData_df,colData(CoCena_Input_HC)[rownames(subsetData_df),"Merged"])
colnames(subsetData_df)[ncol(subsetData_df)] <- "Merged"
mean_per_condition = aggregate(x = subsetData_df[,-ncol(subsetData_df)],
                               by = list(as.factor(subsetData_df[,"Merged"])),
                               FUN = mean)


rownames(mean_per_condition) <- mean_per_condition$Group.1
mean_per_condition$Group.1 <- NULL
colnames(mean_per_condition) <- rowData(CoCena_Input_HC)[rownames(rowData(CoCena_Input_HC)) %in% colnames(mean_per_condition),"SYMBOL"]
z_values <- mean_per_condition
colnames(z_values) <- colnames(mean_per_condition)
nBins_each <- 5
myColor <- c(colorRampPalette(c("#0077b6", "#d9effa"))(nBins_each),
             "white",
             colorRampPalette(c("#ffe8e8","#D62828"))(nBins_each)
)


names(colorTheme) <-c("wt_cdcdcd_HC","ko_cdcdcd_HC",
                      "wt_cdcdhfd_HC","ko_cdcdhfd_HC",
                      "wt_hfdcdcd_HC","ko_hfdcdcd_HC",
                      "wt_hfdhfdhfd_HC","ko_hfdhfdhfd_HC")

annoCol <- list(group = colorTheme)
col_anno = data.frame(group = names(colorTheme))
rownames(col_anno) = col_anno$group

heatmap <- pheatmap(z_values,
                    color = myColor,
                    scale="column",
                    fontsize_number = 15,
                    cellheight = 10,
                    cellwidth = 10,
                    annotation_colors = annoCol,
                    annotation_row  = col_anno,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    angle_col = 90,
                    cluster_rows = F
)

# ORA - Analysis ----
#Universe HC
universe_entrez_HC <- bitr(rownames(res_HC),
                        fromType="ENSEMBL",
                        toType="ENTREZID",
                        OrgDb="org.Mm.eg.db")$ENTREZID
#Universe KC
universe_entrez_KC <- bitr(rownames(res_KC),
                        fromType="ENSEMBL",
                        toType="ENTREZID",
                        OrgDb="org.Mm.eg.db")$ENTREZID

#Universe KC_wt
universe_entrez_KC_wt <- bitr(rownames(res_KC_wt),
                           fromType="ENSEMBL",
                           toType="ENTREZID",
                           OrgDb="org.Mm.eg.db")$ENTREZID

#Universe KC_wt P0
universe_entrez_KC_wt_P0 <- bitr(rownames(res_KC_wt_P0),
                              fromType="ENSEMBL",
                              toType="ENTREZID",
                              OrgDb="org.Mm.eg.db")$ENTREZID

#Universe KC_wt P0
universe_entrez_KC_wt_E14_5 <- bitr(rownames(res_KC_wt_E14_5),
                                 fromType="ENSEMBL",
                                 toType="ENTREZID",
                                 OrgDb="org.Mm.eg.db")$ENTREZID

## ORA ko hfd cd cd vs wt hfd cd cd ----
HC_ORA_UP <- doOra(
  rownames(res_HC[which(res_HC$log2FoldChange>0 & res_HC$padj <0.1),]), # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=4,
  universe_entrez_HC,
  filename = "ORA_DE/HC_UP" # will be ORA_[filename][type].png
  )

HC_ORA_DOWN <- doOra(
  rownames(res_HC[which(res_HC$log2FoldChange<0 & res_HC$padj <0.1),]), # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=4,
  universe_entrez_HC,
  filename = "ORA_DE/HC_DOWN"
)

KC_ORA_UP <- doOra(
  rownames(res_KC[which(res_KC$log2FoldChange>0 & res_KC$padj <0.1),]), # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=6,
  universe_entrez_KC,
  filename = "ORA_DE/KC_UP"
)
KC_ORA_DOWN <- doOra(
  rownames(res_KC[which(res_KC$log2FoldChange<0 & res_KC$padj <0.1),]), # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=5,
  universe_entrez_KC,
  filename = "ORA_DE/KC_DOWN",
  GoFilter = T # Turn on here as many overlapping terms enriched
)

## ORA wt ----
# Identify DE between all Maternal Lean vs all Maternal Obese
maternal_lean <- c("cdcdcd","cdcdhfd")
maternal_obese <- c("hfdcdcd","hfdhfdcd","hfdcdhfd","hfdhfdhfd")

list_of_all_combination_wt <- list()
for(i in maternal_lean){
  for(j in maternal_obese){
    nameList <- paste0(i,"_vs_",j)
    tmp <- results(de_seq_result_KC_WT,
                         contrast =c("Condition",j,i),
                         alpha = 0.1)
    list_of_all_combination_wt[[nameList]] <- tmp
    list_of_all_combination_wt[[paste0(nameList,"_UP")]] <- rownames(tmp[which(tmp$log2FoldChange>0& tmp$padj <0.1),])
    list_of_all_combination_wt[[paste0(nameList,"_DOWN")]] <- rownames(tmp[which(tmp$log2FoldChange<0 & tmp$padj <0.1),])
  }
}

all_up_in_MaternalObese <- unique(unlist(list_of_all_combination_wt[grepl("UP",names(list_of_all_combination_wt))],use.names = F))

all_up_in_MaternalLean <- unique(unlist(list_of_all_combination_wt[grepl("DOWN",names(list_of_all_combination_wt))],use.names = F))

KC_WT_ORA_UP <- doOra(
  all_up_in_MaternalObese, # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=6,
  universe_entrez_KC_wt,
  filename = "ORA_DE/KC_WT_UP"
)
KC_ORA_DOWN <- doOra(
  all_up_in_MaternalLean, # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=5,
  universe_entrez_KC_wt,
  filename = "ORA_DE/KC_WT_DOWN",
  GoFilter = T # Turn on here as many overlapping terms enriched
)


## ORA wt P0 ----
# Identify DE between all Maternal Lean vs all Maternal Obese
diet_difference <- c("CD","HFD")

tmp <- results(de_seq_result_KC_WT_P0,
               contrast =c("Condition","HFD","CD"),
               alpha = 0.1)

all_up_in_MaternalObese <- rownames(tmp[which(tmp$log2FoldChange>2 & tmp$padj <0.1),])
all_up_in_MaternalLean <- rownames(tmp[which(tmp$log2FoldChange<(-2) & tmp$padj <0.1),])

list_of_all_to_Compare <- list(P0_UP = all_up_in_MaternalObese,
                               P0_DOWN = all_up_in_MaternalLean)

KC_WT_P0_ORA_UP <- doOra(
  all_up_in_MaternalObese, # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=6,
  universe_entrez_KC_wt_P0,
  filename = "ORA_DE/KC_WT_P0_UP"
)

KC_WT_P0_ORA_DOWN <- doOra(
  all_up_in_MaternalLean, # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=6,
  universe_entrez_KC_wt_P0,
  filename = "ORA_DE/KC_WT_P0_DOWN"
)

## ORA wt E14.5 ----
# Identify DE between all Maternal Lean vs all Maternal Obese
diet_difference <- c("CD","HFD")

tmp <- results(de_seq_result_KC_WT_E14_5,
               contrast =c("Condition","HFD","CD"),
               alpha = 0.1)


all_up_in_MaternalObese <- rownames(tmp[which(tmp$log2FoldChange>1.2 & tmp$padj < 0.1),])
all_up_in_MaternalLean <- rownames(tmp[which(tmp$log2FoldChange<(-1.2) & tmp$padj < 0.1),])

list_of_all_to_Compare$E14_5_UP =all_up_in_MaternalObese
list_of_all_to_Compare$E14_5_DOWN =all_up_in_MaternalLean

KC_WT_E14_5_ORA_UP <- doOra(
  all_up_in_MaternalObese, # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=6,
  universe_entrez_KC_wt_E14_5,
  filename = "ORA_DE/KC_WT_E14_5_UP"
)

KC_WT_E14_5_ORA_DOWN <- doOra(
  all_up_in_MaternalLean, # ENSEBML
  type=c("KEGG","GO","HALLMARK"),
  levelGOTerms=6,
  universe_entrez_KC_wt_E14_5,
  filename = "ORA_DE/KC_WT_E14_5_DOWN"
)

## Do Overlap analysis all youngDE ----
library(UpSetR)
X <- upset(fromList(list_of_all_to_Compare), order.by = "freq",text.scale = 3)

# Identify the Unique elements
df2 <- data.frame(gene=unique(unlist(list_of_all_to_Compare)))
df1 <- lapply(list_of_all_to_Compare,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "type")

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(type) %>% 
    pull("type") %>% 
    paste0(collapse = "|")
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()


# Intresting due to oppsoing directions DE
rowData(de_seq_result_KC_WT_E14_5)[subset(df_int,subset = int == "E14_5_DOWN|P0_UP")$gene,"SYMBOL"]
rowData(de_seq_result_KC_WT_E14_5)[subset(df_int,subset = int == "E14_5_UP|P0_DOWN")$gene,"SYMBOL"]

#
output_result_list[["ORA_HC"]] <- list(UP = HC_ORA_UP,
                                       DOWN = HC_ORA_DOWN
                                       )  

output_result_list[["ORA_KC"]] <- list(UP = KC_ORA_UP,
                                       DOWN = KC_ORA_DOWN
                                       ) 

output_result_list[["ORA_KC_WT_P0"]] <- list(UP = KC_WT_P0_ORA_UP,
                                       DOWN = KC_WT_P0_ORA_DOWN
) 

output_result_list[["ORA_KC_WT_E14_5"]] <- list(UP = KC_WT_E14_5_ORA_UP,
                                             DOWN = KC_WT_E14_5_ORA_DOWN
) 


# Ligand-Receptor Analysis ----
lr_mouse_table <- read.table("../data/mouse_lr_pair_fromCellTalk.txt",header = T)

# Receptors for sure expressed in "ko_hfdcdcd_HC" or "wt_hfdcdcd_HC" searched for
# for this filtering on subset of initial data

dds_HC_preFilter[,dds_HC_preFilter$Merged %in% c("ko_hfdcdcd_HC","wt_hfdcdcd_HC")]
dds_HC_highThreshold <- preprocessing(dds_HC_preFilter[,dds_HC_preFilter$Merged %in% c("ko_hfdcdcd_HC","wt_hfdcdcd_HC")],
                        11, # >10
                        protCodingOnly=T,
                        removeConstRows=T,
                        filterPerSample=T
)


DE_LigandReceptor <- doLigandReceptorPlot(
  allDEGenes_ligand = res_KC[which(res_KC$padj < 0.1 & abs(res_KC$log2FoldChange)>2),],
  dds_obj_Ligand = dds_KC,
  allPresent_receptor = dds_HC_highThreshold,
  colorVar = "log2FoldChange",
  adjMatrix_LigandReceptor = lr_mouse_table
  )

allLigands <- as.character(unique(lr_mouse_table$ligand_gene_id))
HC_ORA_UP <- doOra(
  unique(DE_LigandReceptor$ligand_ensembl_gene_id[which(DE_LigandReceptor$LFC_noTrim<0)]), # ENSEBML
  type=c("GO"),
  levelGOTerms=6,
  allLigands,
  filename = "ORA_DE/Ligand_Set" # will be ORA_[filename][type].png
)

# CoCena ----
## HC Co-expression analysis ----
# This was done closely following this:
# https://github.com/MarieOestreich/hCoCena/tree/main/showcase
# main and satellite showcase Notebooks
# The Input is generated in here and saved to data folder

CoCena_Input_HC <- DESeq(dds_HC)

saveRDS(CoCena_Input_HC,"../data/CoCena_Input_HC.rds")

## KC Co-expression analysis ----
### Hif1a ----

CoCena_Input_KC <- DESeq(dds_KC)

saveRDS(CoCena_Input_KC,"../data/CoCena_Input_KC.rds")

### WT ----
CoCena_Input_KC_WT <- DESeq(dds_KC_WT)

saveRDS(CoCena_Input_KC_WT,"../data/CoCena_Input_KC_WT.rds")


# Save everything ----
saveRDS(output_result_list, file = "Transcriptomics_results.rds")

setwd("..")

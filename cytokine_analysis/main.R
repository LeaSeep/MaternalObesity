# Ctyokine Analysis ----
## Overview ----
# This script visualizes the cytokine measurements as well as performs stats testing
# to investigate into changes

## Input ----
# (all in data Folder):
#   - mother (obese vs lean) cytokine measurements
#   - offspring (different diet conditions) cytokine measurements
## Output ----
# (output to script's folder):
#   - Cytokine Heatmap for the Mothers
#   - Cytokine Heatmap for the Offspring

# Data Preparation ----
setwd("cytokine_analysis")
cytokine_mother <- as.data.frame(readxl::read_excel("../data/cytokines_serum_mother.xlsx"))
cytokine_offspring <- as.data.frame(readxl::read_excel("../data/cytokines_serum_offspring.xlsx"))

rownames(cytokine_mother) <- cytokine_mother$`Sample ID`
cytokine_mother$`Sample ID` <- NULL
cytokine_mother$Well <- NULL

rownames(cytokine_offspring)<- paste0(cytokine_offspring$`pg/ml`,1:nrow(cytokine_offspring))
cytokine_offspring$`pg/ml`<-NULL


nBins_each <- 15
myColor <- c(
  #colorRampPalette(c("#0077b6", "#d9effa"))(nBins_each),
  #"white",
  colorRampPalette(c("white","#D62828"))(nBins_each)
)


# Mothers analysis ----
colorTheme <- c("#797979","#8fd091")
names(colorTheme) <-c("CD","HFD")
annoCol <- list(group = colorTheme)
col_anno = data.frame(group = rownames(cytokine_mother))
rownames(col_anno) = col_anno$group
col_anno$group <- gsub(" .*","",col_anno$group)

data2Plot <- t(log10(cytokine_mother+1))
# get order based on intensity
data2Plot<-data2Plot[order(apply(cytokine_mother,2,mean),decreasing = T),]

## Heatmap ----
heatmap <- pheatmap(data2Plot,
                    color = myColor,
                    scale="none",
                    breaks = seq(0,4,length.out=15), # no manual adjustment needed if scaling 
                    legend_breaks = seq(0,4,length.out=5),
                    legend_labels = c("10e0", "10e1", "10e2", "10e3", "10e4"),
                    fontsize_number = 15,
                    cellheight = 10,
                    cellwidth = 10,
                    filename = "Heatmap_Mothers.png",
                    annotation_colors = annoCol,
                    annotation_col  = col_anno,
                    cluster_cols = F,
                    cluster_rows = F,
                    angle_col = 90
)
## Testing ----
# test If there is a difference amon the groups
data2test <- cytokine_mother
data2test$group <- gsub(" [0-9]*$","",rownames(data2test))
long <- reshape2::melt(data2test)

pVal_mother <- c()
for(i in unique(long$variable)){
  tmp_df <-subset(long,variable==i)
  tmp=wilcox.test(tmp_df$value ~ tmp_df$group)
  pVal_mother <- c(pVal_mother, setNames(tmp$p.value,i))
}
table(pVal_mother<0.05)
pVal_mother[pVal_mother<0.05]
table(p.adjust(pVal_mother,method = "fdr")<0.1)


# Offspring ----

colorTheme <- c("#c6c6c6","#606060","#c32b38","#fee2d1","#8a4094","#ebdeec")
names(colorTheme) <-c("CDCDCD",
                      "CDCDHFD",
                      "HFDCDCD",
                      "HFDHFDCD",
                      "HFDCDHFD",
                      "HFDHFDHFD")

annoCol <- list(group = colorTheme)
col_anno = data.frame(group = rownames(cytokine_offspring))
rownames(col_anno) = col_anno$group

col_anno$group <- gsub("[0-9].*","",col_anno$group)

data2Plot <- t(log10(cytokine_offspring+1))
# get order based on intensity
data2Plot<-data2Plot[order(apply(cytokine_offspring,2,mean),decreasing = T),]
## Heatmap ----
heatmap <- pheatmap(data2Plot,
                    color = myColor,
                    scale="none",
                    breaks = seq(0,4,length.out=15), # no manual adjustment needed if scaling 
                    legend_breaks = seq(0,4,length.out=5),
                    legend_labels = c("10e0", "10e1", "10e2", "10e3", "10e4"),
                    fontsize_number = 15,
                    cellheight = 10,
                    cellwidth = 10,
                    filename = "Heatmap_Offspring.png",
                    annotation_colors = annoCol,
                    annotation_col  = col_anno,
                    cluster_cols = F,
                    cluster_rows = F,
                    angle_col = 90
)

## Testing ----
# test If there is a difference amon the groups
data2test <- cytokine_offspring
data2test$group <- gsub("[0-9]*$","",rownames(data2test))
long <- reshape2::melt(data2test)

pVal_offspring <- c()
for(i in unique(long$variable)){
  tmp=kruskal.test(formula = value ~ group,data=subset(long,variable==i))
  pVal_offspring <- c(pVal_offspring, setNames(tmp$p.value,i))
}

table(pVal_offspring<0.05)
pVal_offspring[pVal_offspring<0.05]
table(p.adjust(pVal_offspring,method = "fdr")<0.1)


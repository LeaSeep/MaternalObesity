# DecoupleR TF analysis ----
## Overview ----
# Performing TF analysis in the wt data, measured at P0 (birth) to determine
# which regulatory program seems to be causing further developmental issues in
# the offspring born to obese mothers

# Analysis closely follows https://saezlab.github.io/decoupleR/articles/tf_bk.html

## Input ----
# (all in data Folder):
#   - dds Object for KC_WT_P0 - created within main.R within the transcriptome_analysis
#   - DESeq result Object for KC_WT_P0 - created within main.R within the transcriptome_analysis


## Output----
# (output to script's folder):
#   - Score-Plot of TF analysis
#   - Volcano Plot with marked Hif1a targets
#   - Heatmap of DEGs (including Hif1a targets)

# Prep Data ----
## prep libaries
setwd("transcriptome_analysis")

library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)

source("../utils/doComBat.R")


# read in data
dds_data <- readRDS("../data/dds_KC_WT_P0.rds")
res_data <- readRDS("../data/res_KC_WT_P0.rds")

## count data ----
# Perfom batch correction
vst_dds_data <- vst(dds_data,blind=T)
correctedObj <- doBatchCorrection(
  SumExp_obj = vst_dds_data,
  design_factor = "Condition",
  batch_factor = "date_batch"
)

corrected_counts <-  as.data.frame(assay(correctedObj))
rownames_test <- rownames(corrected_counts)
corrected_counts$gene <- rowData(dds_data)[rownames(corrected_counts),"SYMBOL"]
corrected_counts <- corrected_counts[!duplicated(corrected_counts$gene),]
rownames(corrected_counts)<-NULL

## DE result data ----
result_DE <- res_data
result_DE$ID_ens <- rownames(result_DE)
result_DE <- as.data.frame(result_DE)
result_DE$SYMBOL <- rowData(dds_data)[rownames(result_DE),"SYMBOL"]
test_res <- result_DE[!duplicated(result_DE$SYMBOL),]
rownames(test_res) <- NULL

deg <- test_res %>%
  select(SYMBOL, log2FoldChange, stat, pvalue,padj) %>% 
  filter(!is.na(stat)) %>% 
  column_to_rownames(var = "SYMBOL") %>%
  as.matrix()

deg <- as.data.frame(deg)
deg <- deg[!is.na(deg$log2FoldChange),]
deg_all_forVolcano <- deg

deg <- deg[deg[,"padj"]<0.1,]
deg <- deg[abs(deg[,"log2FoldChange"])>2,]
deg <- deg[!is.na(deg$pvalue),]
print(paste0(nrow(deg)," genes are DE under specified criteria"))

# Retrieve Interaction data ----
# takes a while!
interactions <- decoupleR::get_collectri(organism=10090, split_complexes=FALSE)
net <- interactions

# Run ulm ----
contrast_acts <- run_ulm(mat=deg[, 'stat', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor',
                         minsize = 5)


# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))

n_tfs <- 10
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)

f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)


# Plot TF vs their regulatory score ----
score_plot <- ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) +
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred",
                       mid = "whitesmoke", midpoint = 0) +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x =
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold")) +
  xlab("TFs")

ggsave(filename = paste0("KC_WT_P0_TFscorePlot_",Sys.Date(),".png"), plot=score_plot)
ggsave(filename = paste0("KC_WT_P0_TFscorePlot_",Sys.Date(),".svg"), plot=score_plot)

# Display Hif1a targets ----
## Volcano all DE ----
# of all Hif1a targets just label the DE intersection
Hif1aTargets <- unique(net[net$source == "Hif1a","target"])
inter <- rownames(deg)[rownames(deg) %in%Hif1aTargets$target]
deg_all_forVolcano

deg_all_forVolcano$Hif1aDE_Target <- "no"
deg_all_forVolcano[inter,"Hif1aDE_Target"]<-"yes"
deg_all_forVolcano$ID <- NA
deg_all_forVolcano[inter,"ID"] <- inter
deg_all_forVolcano$color <- "grey"
deg_all_forVolcano[rownames(deg),"color"] <- "black"

deg_all_forVolcano[!is.na(deg_all_forVolcano$ID),"color"] <- "darkred"

Volcano <- ggplot(deg_all_forVolcano, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha=0.8,aes(shape=Hif1aDE_Target, color = color)) +
  scale_shape_manual(name="Hif1a Target",labels=c("no","yes"),values=c(20,8))+
  scale_colour_manual(name="color",labels=c("DE","Hif1a Target","non-sig"),values = c("black","darkred","grey")) +
  theme_bw() +
  geom_vline(xintercept = -2, linetype = 'dotted',color="black") +
  geom_vline(xintercept = 2, linetype = 'dotted',color="black") +
  geom_hline(yintercept = -log10(max(deg$pvalue)), linetype = 'dotted',color="black") +
  geom_label_repel(aes(label = ID),color="darkred", size=7)+
  theme(text = element_text(size=20))

ggsave(filename = paste0("KC_WT_P0_VolcanoTargetsMarked_",Sys.Date(),".png"), plot=Volcano)
ggsave(filename = paste0("KC_WT_P0_VolcanoTargetsMarked_",Sys.Date(),".svg"), plot=Volcano)


## Heatmap all DE ----
prepForHeatmap <- corrected_counts[corrected_counts$gene%in%rownames(deg),]
rownames(prepForHeatmap) <- prepForHeatmap$gene
prepForHeatmap$gene <- NULL

palette_length = 100
my_color = colorRampPalette(c("darkblue", "white","red"))(palette_length)
colorTheme <- c("#8fd091","#797979")

names(colorTheme) <- unique(colData(dds_data)[,"Merged",drop=T])
colorTheme <- colorTheme[1:length(unique(colData(dds_data)[,"Merged",drop=T]))]
annoCol <- list(group = colorTheme)

col_anno <- colData(dds_data)[,"Merged",drop=F]
col_anno = as.data.frame(col_anno)
colnames(col_anno) = "group"

pheatmap(prepForHeatmap,
         border_color = NA,
         color = my_color,
         annotation_colors = annoCol,
         annotation_col  = col_anno,
         #breaks = my_breaks,
         cluster_rows = T,
         cluster_cols = T,
         fontsize_number = 15,
         cellheight = 10,
         cellwidth = 10,
         scale = "row",
         filename = "KC_WT_P0_Heatmap_DEG.png")

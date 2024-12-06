# Proteomics ----
# read in proteomics data and put to SumExp object
data_prot <- xlsx::read.xlsx("../data/t-test_HFD-WTvsKO_FDR-0_1.xlsx",sheetName = "LS_dataForHeatmap")

# two duplicates take one with highest rowsum
# Find non-unique values
non_unique_values <- data_prot$NA.[duplicated(data_prot$NA.) | duplicated(data_prot$NA., fromLast = TRUE)]
# Get unique non-unique values
non_unique_values <- unique(non_unique_values)

removeIDX <- c()
removeIDX <- c(removeIDX,names(which.min(rowSums(data_prot[grepl("Gbp4",data_prot$NA.),-1]))))
removeIDX <- c(removeIDX,names(which.min(rowSums(data_prot[grepl("Wdr20",data_prot$NA.),-1]))))

data_prot <- data_prot[-as.numeric(removeIDX),]

# remove all NA rows
data_prot <- data_prot[-which(is.na(data_prot$NA.)),]

rownames(data_prot) <- data_prot$NA.
data_prot$NA. <- NULL

data_sampleAnno <- as.data.frame(readxl::read_excel("../data/t-test_HFD-WTvsKO_FDR-0_1.xlsx",sheet = "LS_sampleAnno"))
rownames(data_sampleAnno) <- data_sampleAnno$...1
data_sampleAnno$...1 <- NULL

data_rowAnno <-  as.data.frame(readxl::read_excel("../data/t-test_HFD-WTvsKO_FDR-0_1.xlsx",sheet = "LS_entitieAnno"))
data_rowAnno <- data_rowAnno[!duplicated(data_rowAnno$...1),]

data_rowAnno <- data_rowAnno[which(data_rowAnno$...1 %in% rownames(data_prot)),]
rownames(data_rowAnno) <- data_rowAnno$...1
data_rowAnno$...1 <- NULL

proteomics <- SummarizedExperiment::SummarizedExperiment(assays = as.matrix(data_prot),
                                                         colData = data_sampleAnno,
                                                         rowData = data_rowAnno)

## get Transport ----
intrestingTerms <- c("lipid_transport","lipoprotein","lipoprotein_metabolic_process","cross")
library(pheatmap)
library(SummarizedExperiment)

i ="cross"
plot <- list()
for(i in intrestingTerms){
  proteomics_tramsport <- proteomics[rowData(proteomics)[,i]==1,]
  data <- proteomics_tramsport
  subsetData_df <- as.data.frame(t(assay(data)))
  
  groupFactor = "merged"
  
  subsetData_df <- cbind(subsetData_df,colData(data)[rownames(subsetData_df),groupFactor])
  
  colnames(subsetData_df)[(ncol(subsetData_df)-length(groupFactor)+1):ncol(subsetData_df)] <- groupFactor
  z_values <- subsetData_df[,1:(ncol(subsetData_df)-length(groupFactor))]
  
  rowData(data)
  colnames(z_values) <- rowData(data)[colnames(z_values),"geneName"]
  nBins_each <- 5
  myColor <- c(colorRampPalette(c("#0077b6", "#d9effa"))(nBins_each),
               "white",
               colorRampPalette(c("#ffe8e8","#D62828"))(nBins_each)
  )
  myColor <- viridis_pal(option = "inferno")(nBins_each*2+1)
  
  colorTheme = c("#a6cee3","#1f78b4",
                 "#fdbf6f","#ff7f00")
  
  names(colorTheme) <- c("WT CD","KO CD","WT HFD","KO HFD")
  
  #names(colorTheme) <- unique(subsetData_df[,groupFactor,drop=T])
  colorTheme <- colorTheme[1:length(unique(subsetData_df[,groupFactor,drop=T]))]
  annoCol <- list(group = colorTheme)
  
  col_anno <- subsetData_df[,groupFactor,drop=F]
  col_anno = as.data.frame(col_anno)
  colnames(col_anno) = "group"
  
  annoCol <- list(group = colorTheme)
  annotation_row_df <- data.frame(row.names=rownames(subsetData_df),group=subsetData_df[,groupFactor])
  
  toplot <- as.data.frame(t(z_values))
  toplot <- as.data.frame(t(apply(toplot,1,scale)))
  colnames(toplot) <- rownames(z_values)
  
  
  
  plot[[i]] <-  pheatmap(toplot,
                         color = myColor,
                         scale = "none",
                         cutree_cols = 2,
                         cutree_rows = 4,
                         # display_numbers = as.data.frame(t(result_sig)), # needs to be a matrix of same dim as LFC Table
                         #fontsize_number = 15,
                         #cellheight = 10,
                         #cellwidth = 10,
                         #filename = givenFilename,
                         show_rownames = F,
                         annotation_colors = annoCol,
                         annotation_col  = annotation_row_df,
                         # treeheight_row = 0,
                         # treeheight_col = 0,
                         # cutree_rows =8,
                         # main = "robust BBS1 & Lean Metabolites",
                         # angle_col = 90,
                         cluster_rows = T,
                         cluster_cols = T,
                         clustering_distance_cols = "correlation",
                         clustering_distance_rows = "correlation",
                         main=i
  )
}

plot$lipid_transport

# take proteomics and subset
allProteinsOfInterest <- xlsx::read.xlsx("../data/venn_diagram_info_HH_selected for Heatmap highlight.xlsx",
                                         sheetName = "Sheet1",header =F)[,c("X4","X5","X6")]
allProteinsOfInteresxt_LookUp <- allProteinsOfInterest[-(1),]
colnames(allProteinsOfInteresxt_LookUp) <- allProteinsOfInteresxt_LookUp[1,]
allProteinsOfInteresxt_LookUp <- allProteinsOfInteresxt_LookUp[-1,]

#subset proteomics
allProts <- unique(unlist(allProteinsOfInteresxt_LookUp,use.names = F))[!is.na(unique(unlist(allProteinsOfInteresxt_LookUp,use.names = F)))]
length(allProts)

proteomics_subset <- proteomics[rownames(proteomics) %in% allProts,]
annotation_row_df <- data.frame(row.names=colnames(proteomics_subset),group=proteomics_subset$merged)

toplot <- as.data.frame((assay(proteomics_subset)))
# z-scale
toplot <- as.data.frame(t(apply(toplot,1,scale)))
colnames(toplot) <- colnames(as.data.frame((assay(proteomics_subset))))

# try to identify clusters
# show mean over samples
toPlot <- as.data.frame(t(toplot))

toPlot$group <- annotation_row_df[rownames(toPlot),"group"]

group_means <- toPlot %>%
  group_by(group) %>%
  summarise_all(mean)



group_means <- as.data.frame(group_means)
rownames(group_means) <- group_means$group
group_means$group <- NULL




annotation_row_df <- data.frame(row.names=rownames(group_means),group=rownames(group_means))

# to label 
toLable <- xlsx::read.xlsx("../data/venn_diagram_info_HH_selected for Heatmap highlight.xlsx",
                           sheetName = "highlight",header = F)[,"X1"]

table(toLable %in% colnames(group_means))
toLable[!toLable %in% colnames(group_means)]


# # Set rows you don't want to label to NA or ""
# colnames(group_means)[!(colnames(group_means) %in% toLable)] <- ""
# 
# 
# # Create a vector of row names, labeling only the rows of interest
# row_labels <- rownames(t(group_means))
# # Set rows you don't want to label to NA or ""
# row_labels[!(rownames(t(group_means)) %in% toLable)] <- ""
# 
# # Create a new annotation for the rows you want to label
# row_border_anno <- ifelse(row_labels == "", "None", "Label")
# 
# # Convert this to a data frame for annotation_row
# row_border_anno <- data.frame(Border = row_border_anno)
# rownames(row_border_anno) <- rownames(group_means)


pheatmap(t(group_means),
         color = myColor,
         scale = "none",
         #cutree_cols = 2,
         cutree_rows = 6,
         # display_numbers = as.data.frame(t(result_sig)), # needs to be a matrix of same dim as LFC Table
         #fontsize_number = 15,
         # cellheight = 5,
         #cellwidth = c(1,2,3,4),
         #filename = givenFilename,
         show_rownames = T,
         annotation_colors = annoCol,
         annotation_col  = annotation_row_df,
         #annotation_row = row_border_anno,
         # treeheight_row = 0,
         # treeheight_col = 0,
         # cutree_rows =8,
         # main = "robust BBS1 & Lean Metabolites",
         # angle_col = 90,
         cluster_rows = T,
         cluster_cols = T,
         # clustering_distance_cols = "correlation",
         # clustering_distance_rows = "correlation",
)

# Step 1: Generate clusters
row_clusters <- cutree(hclust(dist(t(group_means))), k = 6)

# Step 2: Create a color mapping for clusters
cluster_colors <- rainbow(6)  # Generates 6 distinct colors
names(cluster_colors) <- as.character(1:6)

# Step 3: Create an annotation data frame
row_annotation <- data.frame(Cluster = factor(row_clusters))

xlsx::write.xlsx(row_annotation,file = "../results/cluster_annotation_protein_cross.xlsx")

forHao <- row_annotation[order(row_annotation$Cluster),]

annoCol[["Cluster"]] <- cluster_colors


pheatmap(t(group_means),
         color = myColor,
         scale = "none",
         cutree_rows = 6,
         show_rownames = F,
         annotation_colors = annoCol,
         annotation_col  = annotation_row_df,
         annotation_row = row_annotation,
         # annotation_colors = row_annotation_colors, # This adds color to clusters
         cluster_rows = T,
         cluster_cols = T,
         main = i
)

toLable_clusters <- row_clusters[colnames(group_means) %in% toLable]

# for plotting =>
rwocolor_df <- as.data.frame(toLable_clusters)
colnames(rwocolor_df) <- "Cluster"
# Extract the relevant clusters for the genes in `toLable`
unique_clusters <- unique(toLable_clusters)

# Subset the group_means matrix for the genes in 'toLable'
toLable_data <- group_means[,colnames(group_means) %in% toLable ]

# You can also sort these genes based on their clusters if you want
toLable_data <- toLable_data[,order(toLable_clusters)]

# Load necessary package
library(gridExtra)

# Main heatmap
main_heatmap <- pheatmap(t(group_means),
                         color = myColor,
                         scale = "none",
                         cutree_rows = 6,
                         show_rownames = F,
                         annotation_colors = annoCol,
                         annotation_col  = annotation_row_df,
                         annotation_row = row_annotation,
                         cluster_rows = T,
                         cluster_cols = T,
                         main = i
)

# Smaller heatmap for the genes in `toLable`
label_heatmap <- pheatmap(t(toLable_data),
                          color = myColor,
                          scale = "none",
                          cellheight = 10,
                          cellwidth = 20,
                          show_rownames = T,
                          show_colnames = T,
                          cluster_rows = F,
                          cluster_cols = T,
                          annotation_colors = annoCol,
                          annotation_col  = annotation_row_df,
                          annotation_row = rwocolor_df,
                          main = "Highlighted Genes"
)

# Arrange the two heatmaps side by side
grid.arrange(main_heatmap$gtable, label_heatmap$gtable, ncol = 2)



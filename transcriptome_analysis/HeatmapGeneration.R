# Heatmap based on COCENA ----
selectBaseData <- c("KO") # or KO

if (selectBaseData=="WT"){
  CoCena_Input_KC <- readRDS("../data/CoCena_Input_KC_WT.rds") # WT
  colorTheme <- c("#c6c6c6","#606060","#c12c38","#e0775f","#f3b694","#fce2d0")
  names(colorTheme) <- c("cdcdcd",
                         "cdcdhfd",
                         "hfdcdcd",
                         "hfdcdhfd",
                         "hfdhfdcd",
                         "hfdhfdhfd")
}else{
  CoCena_Input_KC <- readRDS("../data/CoCena_Input_KC.rds")
  colorTheme = c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
                 "#fdbf6f","#ff7f00","#fb9a99","#e31a1c")
  names(colorTheme) <- c("wt_cdcdcd_KC","ko_cdcdcd_KC",
                         "wt_cdcdhfd_KC","ko_cdcdhfd_KC",
                         "wt_hfdcdcd_KC","ko_hfdcdcd_KC",
                         "wt_hfdhfdhfd_KC","ko_hfdhfdhfd_KC")
}

vst_data <- counts(CoCena_Input_KC, normalized=TRUE)
data2sum <- vst_data

union_selected <- unique(DE_LigandReceptor$ligand_ensembl_gene_id)
length(union_selected)

subsetData <- data2sum[rownames(data2sum) %in% union_selected,]
nrow(subsetData)
subsetData_df <- as.data.frame(t(subsetData))
subsetData_df <- cbind(subsetData_df,colData(CoCena_Input_KC)[rownames(subsetData_df),"Merged"])
colnames(subsetData_df)[ncol(subsetData_df)] <- "Merged"


mean_per_condition = aggregate(x = subsetData_df[,-ncol(subsetData_df)],
                               by = list(as.factor(subsetData_df[,"Merged"])),
                               FUN = mean)


rownames(mean_per_condition) <- gsub("_Adult","",mean_per_condition$Group.1)
mean_per_condition$Group.1 <- NULL




colnames(mean_per_condition) <- rowData(CoCena_Input_KC)[rownames(rowData(CoCena_Input_KC)) %in% colnames(mean_per_condition),"SYMBOL"]


#z_values <- as.data.frame(t(apply(mean_per_condition,1,scale)))
z_values <- mean_per_condition
colnames(z_values) <- colnames(mean_per_condition)
nBins_each <- 5
myColor <- c(colorRampPalette(c("#0077b6", "#d9effa"))(nBins_each),
             "white",
             colorRampPalette(c("#ffe8e8","#D62828"))(nBins_each)
)


annoCol <- list(group = colorTheme)
col_anno = data.frame(group = names(colorTheme))
rownames(col_anno) = col_anno$group

orderH <- c("Egf","Camp","Bmp2","Cxcl12","Clu","Angptl3","Hrg",
            "F2","F9","F8","Inhbe","Fgg","Ahsg","Fgb","Cp","Apoa1",
            "Vtn","Hc","Apob","Fga","Inhbc","Gc","Lrig1","Col18a1")
orderH[!orderH %in% colnames(z_values)]
orderJ <- orderH[orderH %in% colnames(z_values)]

z_values <- z_values[,orderJ]

#Camp,Sema5a

heatmap <- pheatmap(z_values,
                    color = myColor,
                    scale="column",
                    # breaks = myBreaks, # no manual adjustment needed if scaling 
                    # display_numbers = as.data.frame(t(result_sig)), # needs to be a matrix of same dim as LFC Table
                    fontsize_number = 15,
                    cellheight = 10,
                    cellwidth = 10,
                    #filename = givenFilename,
                    annotation_colors = annoCol,
                    annotation_row  = col_anno,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    main = "DE Ligands Set representatives in the wt data",
                    angle_col = 90,
                    cluster_rows = F,
                    cluster_cols = F
)
# Heatmap based on DE results ----
colorTheme = c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
               "#fdbf6f","#ff7f00","#fb9a99","#e31a1c")
names(colorTheme) <- c("wt_cdcdcd_HC","ko_cdcdcd_HC",
                       "wt_cdcdhfd_HC","ko_cdcdhfd_HC",
                       "wt_hfdcdcd_HC","ko_hfdcdcd_HC",
                       "wt_hfdhfdhfd_HC","ko_hfdhfdhfd_HC")

res_HC <- results(de_seq_result_HC,
                  contrast = c("Merged","wt_hfdcdcd_HC","wt_cdcdcd_HC"), #c("Merged","ko_hfdcdcd_HC","wt_hfdcdcd_HC"),
                  alpha = 0.1)
rowData(de_seq_result_HC)[which(rowData(de_seq_result_HC)[,"SYMBOL"]%in% "Pnp"),"SYMBOL"] <- c("Pnp_1","Pnp_2")

for(k in c("Chemokines","Interleukins","Cytokines")){
  tmp <- paste0("../data/ImmPort_",k,".txt")
  list <- read.table(tmp,header = T)
  nrow(list)
  expressedReceptors <- rowData(de_seq_result_HC)[,"SYMBOL"]
  
  idx_match <- c()
  found_from_list <- c()
  found_in_syn <- c()
  for(i in 1:nrow(list)){
    if(any(grepl(paste0("^",list$Symbol[i],"$"),expressedReceptors,ignore.case = T))){
      idx_match <- c(idx_match,which(grepl(paste0("^",list$Symbol[i],"$"),expressedReceptors,ignore.case = T)))
      found_from_list <- c(found_from_list,i)
    }else{
      tmp <- trimws(unlist(strsplit(list$Synonyms[i],split = ",")))
      if(any(tmp != "-")){
        for(j in tmp){
          if(j !="Pnp"){
            if(any(grepl(paste0("^",j,"$"),expressedReceptors,ignore.case = T))){
              idx_match <- c(idx_match,which(grepl(paste0("^",j,"$"),expressedReceptors,ignore.case = T)))
              found_from_list <- c(found_from_list,i)
              found_in_syn <- c(found_in_syn,j)
            }
          }

        }
      }
    }
  }
  length(unique(idx_match))
  idx_match<-unique(idx_match)
  
  # in chrm
  
  chemokines <- as.data.frame(counts(de_seq_result_HC,normalized = T))[idx_match,]
  rownames(chemokines) <- rowData(de_seq_result_HC)[rownames(chemokines),"SYMBOL"]
  chemokines <- as.data.frame(t(chemokines))
  
  
  subsetData_df <- cbind(chemokines,colData(de_seq_result_HC)[rownames(chemokines),"Merged"])
  colnames(subsetData_df)[ncol(subsetData_df)] <- "Merged"
  
  mean_per_condition = aggregate(x = subsetData_df[,-ncol(subsetData_df)],
                                 by = list(as.factor(subsetData_df[,"Merged"])),
                                 FUN = mean)
  
  rownames(mean_per_condition) <- mean_per_condition$Group.1
  mean_per_condition$Group.1 <- NULL
  
  z_values <- mean_per_condition
  colnames(z_values) <- colnames(mean_per_condition)
  nBins_each <- 5
  myColor <- c(colorRampPalette(c("#0077b6", "#d9effa"))(nBins_each),
               "white",
               colorRampPalette(c("#ffe8e8","#D62828"))(nBins_each)
  )
  
  
  annoCol <- list(group = colorTheme)
  col_anno = data.frame(group = names(colorTheme))
  rownames(col_anno) = col_anno$group
  
  heatmap <- pheatmap(z_values,
                      color = myColor,
                      scale="column",
                      # breaks = myBreaks, # no manual adjustment needed if scaling 
                      # display_numbers = as.data.frame(t(result_sig)), # needs to be a matrix of same dim as LFC Table
                      fontsize_number = 15,
                      cellheight = 10,
                      cellwidth = 10,
                      #filename = givenFilename,
                      annotation_colors = annoCol,
                      annotation_row  = col_anno,
                      treeheight_row = 0,
                      treeheight_col = 0,
                      main = k,
                      angle_col = 90,
                      cluster_rows = T,
                      cluster_cols = T
  )
  
}

# Heatmap DE Ligand ----
unique(DE_LigandReceptor$ligand_ensembl_gene_id[which(DE_LigandReceptor$LFC_noTrim<0)])

# KC ko
colorTheme = c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
               "#fdbf6f","#ff7f00","#fb9a99","#e31a1c")
names(colorTheme) <- c("wt_cdcdcd_KC","ko_cdcdcd_KC",
                       "wt_cdcdhfd_KC","ko_cdcdhfd_KC",
                       "wt_hfdcdcd_KC","ko_hfdcdcd_KC",
                       "wt_hfdhfdhfd_KC","ko_hfdhfdhfd_KC")

LigandSet <- as.data.frame(counts(de_seq_result_KC,normalized =T)[unique(DE_LigandReceptor$ligand_ensembl_gene_id[which(DE_LigandReceptor$LFC_noTrim<0)]),])

rownames(LigandSet) <- rowData(de_seq_result_KC)[rownames(LigandSet),"SYMBOL"]
LigandSet <- as.data.frame(t(LigandSet))


subsetData_df <- cbind(LigandSet,colData(de_seq_result_KC)[rownames(LigandSet),"Merged"])
colnames(subsetData_df)[ncol(subsetData_df)] <- "Merged"

mean_per_condition = aggregate(x = subsetData_df[,-ncol(subsetData_df)],
                               by = list(as.factor(subsetData_df[,"Merged"])),
                               FUN = mean)

rownames(mean_per_condition) <- mean_per_condition$Group.1
mean_per_condition$Group.1 <- NULL

z_values <- mean_per_condition
colnames(z_values) <- colnames(mean_per_condition)
nBins_each <- 5
myColor <- c(colorRampPalette(c("#0077b6", "#d9effa"))(nBins_each),
             "white",
             colorRampPalette(c("#ffe8e8","#D62828"))(nBins_each)
)


annoCol <- list(group = colorTheme)
col_anno = data.frame(group = names(colorTheme))
rownames(col_anno) = col_anno$group

heatmap <- pheatmap(z_values,
                    color = myColor,
                    scale="column",
                    # breaks = myBreaks, # no manual adjustment needed if scaling 
                    # display_numbers = as.data.frame(t(result_sig)), # needs to be a matrix of same dim as LFC Table
                    fontsize_number = 15,
                    cellheight = 10,
                    cellwidth = 10,
                    #filename = givenFilename,
                    annotation_colors = annoCol,
                    annotation_row  = col_anno,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    main = "DE Ligand Set over all Conditions",
                    angle_col = 90,
                    cluster_rows = T,
                    cluster_cols = T
)

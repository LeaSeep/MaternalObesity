# additional visualisations

## Pre integration ----
library(ggplot2)
DimPlot(wt_scMulti_CD, 
        reduction = "umap.rna", 
        label = TRUE,
        label.size = 2.5,
        repel = TRUE) + ggtitle(paste0("RNA - log and scaled - ", unique(wt_scMulti_CD$orig.ident)))


DefaultAssay(wt_scMulti_CD) <- "RNA"
FeaturePlot(wt_scMulti_CD, features = "nCount_RNA")
FeaturePlot(wt_scMulti_CD, features = "nFeature_RNA")


DimPlot(wt_scMulti_HFD, 
        reduction = "umap.rna", 
        label = TRUE,
        label.size = 2.5,
        repel = TRUE) + 
  ggtitle(paste0("RNA - log and scaled - ", unique(wt_scMulti_HFD$orig.ident)))
DefaultAssay(wt_scMulti_HFD) <- "RNA"
FeaturePlot(wt_scMulti_HFD, features = "nCount_RNA")
FeaturePlot(wt_scMulti_HFD, features = "nFeature_RNA")


# Indetifying cluster annotation
DefaultAssay(int.allCond) <- "RNA"

svglite::svglite("DotPlot_MarkerOverCellTypes.svg",width = 15, height = 10)
DotPlot(int.allCond_woDoublets,
        scale = F,
        features = c(
          "Top2a",
          "Mki67",
          "Cenpp",
          "Cdkn2c",
          
          "Mrc1",
          "Pecam1",
          "Cd38",
          "Ptprb",
          "F8",
          
          
          "Ebf1",
          "Ighd",
          
          "Spp1",
          "Epcam",
          "Cftr",
          "Sox9",
          "Ppp2r2b",
          "Bmp5",
          "Nrxn1",
          "Reln",
          "Col14a1",
          "Pdgfrb",
          
          "Cd74",
          "Btla",
          "Irf8",
          "H2-Eb1",
          "H2-Aa",
          "Flt3",
          
          "Skap1",
          "Cd226",
          "Il2rb",
          "Trbc2",
          "Bcl11b",
          
          "Adgre1",
          "Vsig4",
          "Clec4f",
          "Timd4",
          "Csf1r",
          
          "Alb",
          "Arg1",
          "Pck1",
          "Hc"
          # "Ccr2",
          # "Flt3",
          # "H2-Eb1"
        ),dot.min = 0.1
        
)+RotatedAxis()+
  theme(aspect.ratio = 0.5)
dev.off()


DimPlot(int.allCond_myeloid, 
        reduction = "umap.myeloid.wnn", 
        label = TRUE, label.size = 5, repel = TRUE,
        split.by  = "orig.ident") + 
  ggtitle("integrated WNN")
Idents(int.allCond_myeloid) <- int.allCond_myeloid$seurat_clusters
DimPlot(int.allCond_myeloid, 
        reduction = "umap.myeloid.wnn", 
        label = TRUE, 
        label.size = 5,
        split.by="orig.ident",
        repel = TRUE,
        cols = colorSubset) + 
  ggtitle("integrated WNN")


DefaultAssay(int.allCond_myeloid) <- "RNA"
FeaturePlot(int.allCond_myeloid,
            features = c("Apoe"),
            reduction = 'umap.myeloid.wnn',
            split.by = "orig.ident")

DefaultAssay(int.allCond_myeloid) <- "SCT"
DotPlot(int.allCond_myeloid,
        features = c("Apoe"),
        #reduction = 'umap.myeloid.wnn',
        split.by = "orig.ident")


# new colors
colorSubset<-c("#cb997e","#2a9d8f",
               "#8ab17d","#e9c46a",
               "#f4a261","#e76f51")
names(colorSubset) <- paste0("KC_",0:5)

colorAll <- scales::hue_pal()(9)
names(colorAll) <- levels(Idents(int.allCond_woDoublets))


names(colorSubset) <- 0:5
Idents(int.allCond_myeloid_woCluster5) <- int.allCond_myeloid_woCluster5$seurat_clusters
DefaultAssay(int.allCond_myeloid_woCluster5) <- "integrated"
Seurat::DimPlot(int.allCond_myeloid_woCluster5, 
                reduction = "umap.myeloid.wnn", 
                label = TRUE, label.size = 5,
                split.by  = "orig.ident", 
                cols = colorSubset,
                repel = TRUE) + 
  ggtitle("integrated WNN")+
  theme(aspect.ratio = 1)


DefaultAssay(int.allCond_myeloid)<-"SCT"
DotPlot(
  int.allCond_myeloid_woCluster5,
  #  int.allCond_myeloid,
  #split.by = "orig.ident",
  features = c("Adgre1",
               "Timd4",
               "Clec4f",
               "C1qa",
               "H2-Aa",
               "H2-Ab1",
               "Cd74",
               "Pparg",
               "Lilr4b",
               "Lgals3",
               "Sema6d",
               "C6",
               "F8",
               "Adam23",
               "Marco",
               "Il10",
               "Tgfb1",
               "Cx3cr1")
)+RotatedAxis()+ theme(aspect.ratio = 1)


# To get all pVals
list_markers_betweenCond <- list()
int.allCond_myeloid_woCluster5_SCTprepped$merged_cond <- paste0(int.allCond_myeloid_woCluster5_SCTprepped$orig.ident,
                                                     "_",
                                                     int.allCond_myeloid_woCluster5_SCTprepped$seurat_clusters)
table(int.allCond_myeloid_woCluster5_SCTprepped$merged_cond)
Idents(int.allCond_myeloid_woCluster5_SCTprepped) <- int.allCond_myeloid_woCluster5_SCTprepped$merged_cond
for(i in 0:4){
  DE <- FindMarkers(int.allCond_myeloid_woCluster5_SCTprepped,
                    ident.1 = paste0("wtHFDCDCD_",i),
                    ident.2 = paste0("wtCDCDCD_",i),
                    assay = "SCT",
                    min.pct = 0.10,
                    logfc.threshold = 0,
                    return.thresh = 0.01,
                    only.pos =F) 
  list_markers_betweenCond[[paste0("cluster",i)]] <- DE
}


DefaultAssay(int.allCond_myeloid_woCluster5_SCTprepped) <- "RNA"
chosenGenes <- "Apoe"
genesOfInterest_basedOnPeak_2 <- c("Gpnmb","Pck1","Apoe","Apoa1")

for(i in genesOfInterest_basedOnPeak_2){
  chosenGenes <- i
  
  vln_plot <- tryCatch({
    VlnPlot(int.allCond_myeloid_woCluster5, 
            features = chosenGenes, 
            group.by = "seurat_clusters",
            split.by = "orig.ident",
            pt.size = 0.01, 
            combine = T,
            add.noise = T)
  }, error = function(e) {
    message(paste("Error with gene:", chosenGenes, "- skipping this iteration"))
    NULL
  })
  
  if (!is.null(vln_plot)) {
    # Extract the data used for plotting
    plot_data <- vln_plot$data
    plot_data$cluster <- gsub("^.*_","",plot_data$ident)
    # Create the violin plot with statistical testing
    
    my_comparisons <- list( c("wtCDCDCD", "wtHFDCDCD"))
    colorVector <- c(rep(c("#c6c6c6","#c12c38")))
    names(colorVector) <- c("wtCDCDCD", "wtHFDCDCD")
    plot_data$split <- factor(plot_data$split, levels = c("wtCDCDCD", "wtHFDCDCD"))
    
   # take p values from overall search
    p_realDeal <- rlist::list.rbind(lapply(list_markers_betweenCond, function(x) x[grepl(chosenGenes,rownames(x)),]))
    expected_clusters <- paste0("cluster", 0:4)
    
    
    
    if(nrow(p_realDeal)==0){
      print(paste0(chosenGenes," nothing DE"))
      png(paste0("expr_boxplot_",chosenGenes,".png"))
      plot(vln_plot)
      dev.off()
    }else{
      print(paste0(chosenGenes,"DE"))
      # Check which clusters are missing
      missing_clusters <- setdiff(expected_clusters, rownames(p_realDeal))
      
      # Create a data frame for the missing clusters with all columns set to 1
      if (length(missing_clusters) > 0) {
        missing_df <- data.frame(
          p_val = rep(1,length(missing_clusters)),
          avg_log2FC = rep(1,length(missing_clusters)),
          pct.1 = rep(1,length(missing_clusters)),
          pct.2 = rep(1,length(missing_clusters)),
          p_val_adj = rep(1,length(missing_clusters)),
          row.names = missing_clusters
        )
        
        # Combine the original data with the missing clusters
        p_realDeal <- rbind(p_realDeal, missing_df)
      }
      # Sort the data frame by row names (optional)
      p_realDeal <- p_realDeal[order(rownames(p_realDeal)), ]
      
      stat.test <- data.frame(
        split = as.character(0:4),
        cluster= as.character(0:4),
        Gene = chosenGenes,
        .y. = chosenGenes,
        group1 = rep("wtCDCDCD",5), 
        group2 = rep("wtHFDCDCD",5), 
        #group1 = paste0(0:4,"wtCDCDCD"), 
        #group2 = paste0(0:4,"wtHFDCDCD"), 
        p_value = format(p_realDeal$p_val,scientific = T,digits = 2),
        y.position = rep(max(plot_data[,chosenGenes])+0.2,5),
        xmin = c(0.8),
        xmax = c(1.2)
      )
      
      # add x axis position for points
      plot_data$x1.position <- 1.2
      plot_data[plot_data$split=="wtCDCDCD","x1.position"] <- 0.8
      final_plot <- ggplot(plot_data, aes(x = cluster, y = plot_data[,chosenGenes], fill = split)) +
        geom_violin(trim = T,scale = "width") + 
        scale_fill_manual(values = colorVector) +
        #geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
        # geom_text(aes(label = p_value), 
        #           position = position_dodge(0.9), 
        #           vjust = -1.5) +
        theme(legend.position = "none") +
        ylab(paste0(chosenGenes," Expression")) + 
        xlab("Condition")+
        theme_bw() +
        facet_wrap(~cluster,ncol=5, scales = "free_x")+
        stat_pvalue_manual(
          tibble::as_tibble(stat.test),
          label = "p_value",
          xmin = "xmin",  # Specify xmin
          xmax = "xmax",  # You might need to define xmax explicitly based on your groups
          y.position = "y.position",  # Ensure y.position is correctly mapped
          group = "split",  # Explicitly set the group aesthetic
          tip.length = 0.01,  # Adjust as necessary
          step.increase = 0
        ) +
        geom_jitter(aes(x = x1.position), 
                    width = 0.18, 
                    size = 0.5) +
        theme(legend.position = "none",axis.text.x=element_blank(), 
              axis.ticks.x=element_blank(),
              aspect.ratio = 2.2)
      
      stat.test$xmin <- stat.test$xmin+c(0,1,2,3,4)
      stat.test$xmax <- stat.test$xmax+c(0,1,2,3,4)
      stat.test
      
      final_plot <- VlnPlot(int.allCond_myeloid_woCluster5, 
                            features = chosenGenes, 
                            group.by = "seurat_clusters",
                            split.by = "orig.ident",
                            pt.size = 1, 
                            combine = T,
                            add.noise = T)+
        scale_fill_manual(values = colorVector) +
        theme(legend.position = "bottom",
              legend.justification = "center",
              axis.text.x=element_text(angle=0)) +
        stat_pvalue_manual(
          tibble::as_tibble(stat.test),
          label = "p_value",
          xmin = "xmin",  # Specify xmin
          xmax = "xmax",  # You might need to define xmax explicitly based on your groups
          y.position = "y.position",  # Ensure y.position is correctly mapped
          group = "split",  # Explicitly set the group aesthetic
          tip.length = 0.0,  # Adjust as necessary
          step.increase = 0
        ) +ylim(0, max(plot_data[,chosenGenes] +0.25))
      
      
      
      svglite::svglite(paste0("chosen_expr_boxplot_",chosenGenes,".svg"))
      print(final_plot)
      dev.off()
    }
  }
}


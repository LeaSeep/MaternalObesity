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

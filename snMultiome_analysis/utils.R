# utils

RNA_integration <- function(seurat.list){
  features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
  seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", 
                                           anchor.features = features)
  seurat.int <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT")
  seurat.int <- RunPCA(seurat.int, verbose = FALSE)
  seurat.int <- RunUMAP(seurat.int, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_'#,
                        #umap.method = 'umap-learn',
                        #metric = 'correlation'
  )
  
  return(seurat.int)
}


RunATACharmony <- function(seurat){
  seurat <- harmony::RunHarmony(object=seurat,
                                group.by.vars='orig.ident',
                                reduction.use='lsi',
                                assay.use='int_peaks',
                                project.dim=F)
  seurat <- RunUMAP(seurat,dims=2:50,reduction='harmony',reduction.name = "umap.atac",reduction.key = "atacUMAP_")
  seurat <- FindNeighbors(seurat, dims = 1:50,reduction='harmony')
  seurat <- FindClusters(seurat, resolution = 2,reduction='harmony')
  seurat <- BuildClusterTree(object = seurat,reorder = T,reorder.numeric = T,verbose = T)
  
  return(seurat)
}

RunWnnRnaAtac <- function(seurat){
  seurat <- FindMultiModalNeighbors(seurat,
                                    reduction.list = list("pca", "harmony"),
                                    dims.list = list(1:50, 2:50),
                                    k.nn = 50)
  seurat <- RunUMAP(seurat, nn.name = "weighted.nn", reduction.name = "INTwnn.umap", reduction.key = "INTwnnUMAP_", assay = "RNA")
  seurat <- FindClusters(seurat, graph.name = "wsnn", 
                         algorithm = 3, verbose = FALSE,
                         resolution = 0.5)
}

GetLSI <- function(seurat){
  seurat <- FindTopFeatures(seurat, min.cutoff = 10)
  seurat <- RunTFIDF(seurat)
  seurat <- RunSVD(seurat)
  seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 2:50)
  
  return(seurat)
}

FindMarkers2Excel <- function(seurat,
                              idents = "orig.ident",
                              grouping.var = "seurat_clusters",
                              filename="./ConservedMarkerList_integratedWNN.xlsx"){
  #Idents(seurat) <- seurat[,idents]
  allMarkers <- FindAllMarkers(seurat,
                               #ident.1= 1,
                               #ident.2= unique(int.allCond[[idents]])[2,1],
                               #grouping.var = idents,
                               min.pct = 0.10,
                               min.diff.pct=0.20,
                               logfc.threshold = 0.20,
                               return.thresh = 0.01,
                               only.pos = T)
  print("Markers Overview:")
  print(table(allMarkers$cluster))
  ### Create list with markers
  totalNrClusters<-max(as.numeric(names(table(allMarkers$cluster))))
  totalNrClustersPlusOne<-totalNrClusters+1
  markersList <- list()
  
  for(i in 1:totalNrClustersPlusOne){
    clusterNr <- i-1
    tmp <- allMarkers[allMarkers$cluster==clusterNr,]
    tmp$score <- tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
    
    markersList[[i]] <- tmp[order(tmp$score, decreasing=TRUE),]
  }
  names(markersList) <- paste0("cluster",0:totalNrClusters)
  
  ### Write to Excel
  library('openxlsx')
  write.xlsx(markersList, file =filename)
  return(markersList)
}


# Create function to get conserved markers for any given cluster
get_conserved <- function(seurat_clusters){
  FindConservedMarkers(int.allCond,
                       ident.1 = seurat_clusters,
                       grouping.var = "orig.ident", # sample
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = seurat_clusters, .)
}



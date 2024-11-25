# utils
library(Signac)
library(EnsDb.Mmusculus.v79)
library(reticulate)
library('openxlsx')
# takes 10 Input File and returns Seurat object with 
# filtered RNA take counts taken diretly from Seurat
# and new generated Macs call
MacPeakCalling <- function(
    absolut_path = "/Volumes/My_Book/Seep_Lea/current/Lipid+Trans/scRNA_ATAC_data/",
    env_path ="./",
    sample_typ = "wtHFDCDCD",
    QC_plots = T
){
  
  # the 10x hdf5 file contains both data types. 
  inputdata.10x <- Read10X_h5(paste0(absolut_path,sample_typ,"/filtered_feature_bc_matrix.h5"))
  rna_counts <- inputdata.10x$`Gene Expression`
  pbmc <- CreateSeuratObject(counts = rna_counts, assay="RNA")
  
  frag.file <- paste0(absolut_path,sample_typ,"/atac_fragments.tsv.gz")

  annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  
  
  
  # create ATAC assay and add it to the object
  atac_peaks <- inputdata.10x$Peaks
  grange.counts <- StringToGRanges(rownames(atac_peaks), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_peaks <- atac_peaks[as.vector(grange.use), ]
  
  seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  genome(annotation) <- "mm10"
  
  pbmc[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_peaks,
    sep = c(":", "-"),
    fragments = frag.file,
    annotation = annotation
  )
  
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, 
                                               pattern = "^mt-",
                                               assay = "RNA")
  
  DefaultAssay(pbmc) <- "ATAC"
  
  pbmc <- NucleosomeSignal(pbmc)
  pbmc <- TSSEnrichment(pbmc)
  
  
  # filter out low quality cells (this currently filters on the RNA (from before moving to MACS to decrease run time a bit)
  pbmc_filtered <- subset(
    x = pbmc,
    subset = 
  #    nCount_ATAC > 1000 &
      nCount_RNA > 500 & #https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
  #    nucleosome_signal < 2 &
  #    TSS.enrichment > 1 &
      percent.mt < 1 &
      nFeature_RNA > 200 & 
      nFeature_RNA < 2500
    
  )
  ## Filtered usage
  # needs a previously create venv!
  tryCatch(reticulate::use_virtualenv("/root/.virtualenvs/r-reticulate"),
           error = function(e){
             warning(paste0("You need to create venv first. Must be in current dir: ",getwd()))
             return(NULL)
           })
  
  reticulate::py_config()
  DefaultAssay(pbmc_filtered) <- "ATAC"
  peaks <- Signac::CallPeaks(pbmc_filtered,macs2.path ="/root/.virtualenvs/r-reticulate/bin/macs3")
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
  
  
  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(pbmc_filtered),
    features = peaks,
    cells = colnames(pbmc_filtered)
  )
  
  # create a new assay using the MACS3 peak set and add it to the Seurat object
  pbmc_filtered[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = frag.file,
    annotation = annotation
  )
  
  
  # compute nucleosome signal score per cell
  DefaultAssay(pbmc_filtered) <- "peaks"
  
  pbmc_filtered <- NucleosomeSignal(object = pbmc_filtered)
  pbmc_filtered <- TSSEnrichment(object = pbmc_filtered, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  total_fragements <- CountFragments(frag.file)
  rownames(total_fragements) <- total_fragements$CB
  pbmc_filtered$fragments <- total_fragements[colnames(pbmc_filtered), "frequency_count"]
  
  pbmc_filtered <- FRiP(
    object = pbmc_filtered,
    assay = 'peaks',
    total.fragments = 'fragments'
  )
  
  pbmc_filtered$blacklist_fraction <- FractionCountsInRegion(
    object = pbmc_filtered, 
    assay = 'peaks',
    regions = blacklist_mm10
  )
  
  
  pbmc_filtered$blacklist_ratio <- pbmc_filtered$blacklist_fraction
  pbmc_filtered$pct_reads_in_peaks <- pbmc_filtered$FRiP
  pbmc_filtered$orig.ident <- sample_typ
  
  # complete filter
  pbmc_filtered <- subset(
    x = pbmc_filtered,
    subset = 
      nCount_ATAC > 1000 &
      nCount_RNA > 500 & #https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
      nucleosome_signal < 2 &
      TSS.enrichment > 1 &
      percent.mt < 1 &
      nFeature_RNA > 200 & 
      nFeature_RNA < 2500
    
  )
  
  if(QC_plots){
    QC_plots <- list()
    # Visualize QC metrics as a violin plot
    
    QC_plots[["all_Feature_Counts_mt_Vln"]] <- VlnPlot(pbmc, 
                                                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                                                       ncol = 3)
    QC_plots[["filtered_Feature_Counts_mt_Vln"]] <- VlnPlot(pbmc_filtered, 
                                                            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                                                            ncol = 3)
    QC_plots[["filtered_rnaC_mt_Scatter"]] <- FeatureScatter(pbmc_filtered,
                                                             feature1 = "nCount_RNA",
                                                             feature2 = "percent.mt")
    QC_plots[["filtered_rnaC_feature_Scatter"]] <- FeatureScatter(pbmc_filtered,
                                                                  feature1 = "nCount_RNA",
                                                                  feature2 = "nFeature_RNA")
    
    
    QC_plots[["all_allQC"]] <- VlnPlot(object = pbmc,
                                       features = c("nCount_RNA",'nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
                                       pt.size = 0,
                                       ncol = 4)
    
    
    QC_plots[["filtered_allQC"]] <- VlnPlot(object = pbmc_filtered,
                                            features = c('nCount_ATAC', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
                                            pt.size = 0,
                                            ncol = 5)
    
    return(list(seuratObj = pbmc_filtered, QC_Plots = QC_plots))
    
  }else{
    return(list(seuratObj = pbmc_filtered, QC_Plots = NULL))
  }
  
}


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



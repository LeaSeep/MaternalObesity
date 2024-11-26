# Main script for sn-transcriptomcis
# snSeq and ATAC seq Analysis

library(Seurat)
library(Signac)
source("utils.R")
library(AnnotationHub)
library(ge)
library(ggplot2)
library(ggsci)
library(hdf5r)
library(glmGamPoi)
library(CellChat)
library(patchwork)
library(presto)
library(decoupleR)
library(OmnipathR)
library(tibble)
library(tidyr)
library('openxlsx')
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
# SetUp ----
wtCDCDCD_orig <- MacPeakCalling(sample_typ = "wtCDCDCD",
                                absolut_path = "/home/rstudio/program/data/",
                                env_path = "/home/rstudio/program/revision_1/MaternalObesity/snMultiome_analysis/")
wtHFDCDCD_orig <- MacPeakCalling(sample_typ = "wtHFDCDCD",
                                 absolut_path = "/home/rstudio/program/data/")


wt_scMulti_CD <- wtCDCDCD_orig$seuratObj
wt_scMulti_HFD <- wtHFDCDCD_orig$seuratObj

# Integration ----
#(following: https://github.com/rebeccaorourke-cu/Sagerstrom_zebrafish_hindbrain/blob/main/workspace/notebooks/03a_Integrate_3WT_neural_subsets.Rmd)

## RNA Prep ----
# gene expression data processing 
### wtCDCDCD ----
DefaultAssay(wt_scMulti_CD) <- "RNA"
wt_scMulti_CD <- SCTransform(wt_scMulti_CD, vars.to.regress = "percent.mt", verbose = FALSE)
wt_scMulti_CD <- RunPCA(wt_scMulti_CD, verbose = FALSE,
                        reduction.name = 'pca.sct',
                        reduction.key = 'sctPC')
wt_scMulti_CD <- RunUMAP(wt_scMulti_CD, dims = 1:50,
                         reduction = 'pca.sct',
                         reduction.name = 'umap.sct', 
                         reduction.key = 'sctUMAP_')
wt_scMulti_CD[["ATAC"]] <- NULL

### wtHFDCDCD ----
DefaultAssay(wt_scMulti_HFD) <- "RNA"
wt_scMulti_HFD <- SCTransform(wt_scMulti_HFD, vars.to.regress = "percent.mt", verbose = FALSE)
wt_scMulti_HFD <- RunPCA(wt_scMulti_HFD, verbose = FALSE,
                         reduction.name ='pca.sct',
                         reduction.key = 'sctPC')
wt_scMulti_HFD <- RunUMAP(wt_scMulti_HFD, dims = 1:50, 
                          verbose = FALSE,
                          reduction = 'pca.sct',
                          reduction.name = 'umap.sct', 
                          reduction.key = 'sctUMAP_')

wt_scMulti_HFD[["ATAC"]] <- NULL

## Across modality integration ----
combined.peaks <- reduce(x = c(wt_scMulti_CD@assays$peaks@ranges,wt_scMulti_HFD@assays$peaks@ranges))

## Call peaks on combined features----
# peaks => MAC called peaks
# int_peaks => peaks called over combined peak ranges
DefaultAssay(wt_scMulti_CD) <- "peaks"
wt_scMulti_CD.counts <- FeatureMatrix(
  fragments = Fragments(wt_scMulti_CD),
  features = combined.peaks,
  cells = colnames(wt_scMulti_CD)
)
DefaultAssay(wt_scMulti_HFD) <- "peaks"
wt_scMulti_HFD.counts <- FeatureMatrix(
  fragments = Fragments(wt_scMulti_HFD),
  features = combined.peaks,
  cells = colnames(wt_scMulti_HFD)
)

ah = AnnotationHub()
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, verbose = T)

seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


wt_scMulti_CD[["int_peaks"]] <- CreateChromatinAssay(
  counts = wt_scMulti_CD.counts,
  fragments = Fragments(wt_scMulti_CD),
  annotation = annotation#,
  # genome = 'mm10'
)


wt_scMulti_HFD[["int_peaks"]] <- CreateChromatinAssay(
  counts = wt_scMulti_HFD.counts,
  fragments = Fragments(wt_scMulti_HFD),
  annotation = annotation#,
  # genome = 'mm10'
)

#saveRDS(wt_scMulti_CD,"../data/CD_scMUlti_Feb.RDS")
#saveRDS(wt_scMulti_HFD,"../data/HFD_scMulti_Feb.RDS")

## find top features for each ----
# peaks
DefaultAssay(wt_scMulti_CD) <- "int_peaks"
# wt_scMulti_CD <- FindTopFeatures(wt_scMulti_CD, min.cutoff = 10)
# wt_scMulti_CD <- RunTFIDF(wt_scMulti_CD)
# wt_scMulti_CD <- RunSVD(wt_scMulti_CD)

DefaultAssay(wt_scMulti_HFD) <- "int_peaks"
# wt_scMulti_HFD <- FindTopFeatures(wt_scMulti_HFD, min.cutoff = 10)
# wt_scMulti_HFD <- RunTFIDF(wt_scMulti_HFD)
# wt_scMulti_HFD <- RunSVD(wt_scMulti_HFD)

# Up till now:
# This integration strategy will be to use the ATAC-seq peaks from all multiomic samples 
# to create a peak set that will be used to add a Chromatin assay called "int_peaks" to each sample. 

## RNA Integration ----
# Now:
# The multiomic sample can then be RNA-seq integrated which will merge the ATAC-seq assays. 
DefaultAssay(wt_scMulti_CD) <- "SCT"
DefaultAssay(wt_scMulti_HFD) <- "SCT"

list.allCond <- list(wt_scMulti_CD,wt_scMulti_HFD)
rm(wt_scMulti_CD)
rm(wt_scMulti_HFD)
int.allCond <- RNA_integration(list.allCond)

#saveRDS(int.allCond,"../data/int.allCond_Feb.RDS")
int.allCond <- RunUMAP(int.allCond, dims = 1:50, 
                       reduction.name = 'umap.rna.int', 
                       reduction.key = 'rnaintUMAP_')

#saveRDS(int.allCond,"../data/int.allCond_Feb.RDS")

## Harmony integration of ATAC-seq ----
# Will then use harmony to integrate the ATAC-seq based on the "int_peaks" assays 


DefaultAssay(int.allCond) <- "int_peaks" # simple merge of both
int.allCond <- GetLSI(int.allCond)
int.allCond <- RunATACharmony(int.allCond)


## Weighted Nearest Neighbor integration of RNAseq and ATACseq ----
# finally will integrate the Seurat integrated RNA-seq and harmony integrated 
# ATAC-seq with weighted nearest neighbor for final cluster determination.
DefaultAssay(int.allCond) <- "integrated"
set.seed(123)
int.allCond <- RunWnnRnaAtac(int.allCond)

saveRDS(int.allCond,"../data/int.allCond_Feb.RDS")
int.allCond <- readRDS("./int.allCond_Feb.RDS")

int.allCond.rna <- DimPlot(int.allCond, reduction = "umap.rna", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("integrated RNA")
int.allCond.atac <- DimPlot(int.allCond, reduction = "umap.atac", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("integrated ATAC")
int.allCond.wnn <- DimPlot(int.allCond, reduction = "INTwnn.umap", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("integrated WNN")
int.allCond.rna + int.allCond.atac + int.allCond.wnn & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

int.allCond.wnn_split <- DimPlot(int.allCond, 
                                 reduction = "INTwnn.umap", 
                                 label = TRUE, label.size = 5, repel = TRUE,
                                 split.by  = "orig.ident") + ggtitle("integrated WNN")


# Cluster Identification ----
## for integrated ----
DefaultAssay(int.allCond) <- "RNA"
int.allCond$merged <- paste0(int.allCond$orig.ident,"_",int.allCond$seurat_clusters)
MarkerList <- FindMarkers2Excel(
  int.allCond,
  idents = "orig.ident",
  grouping.var = "seurat_clusters",
  filename = "./MarkerList_integratedWNN_final.xlsx")


## Annotation for final UMAP all ----
# NOTE THAT THIS ONLY WORKS  WITH THE PROVIDED SET AS THE CLUSTER NUMBER
# ASSIGNMENT IS NOT FIXED, CHECK THE MARKERS BEFORE ASSIGNING
int.allCond$cellType <- as.character(int.allCond$seurat_clusters)
int.allCond$cellType[int.allCond$cellType %in% c("0","7","16")] <- "LSEC"
int.allCond$cellType[int.allCond$cellType %in% c("1","4","5","8",
                                                 "11","12","14")] <- "hepatocytes"
int.allCond$cellType[int.allCond$cellType %in% c("2","6")] <- "KC"
int.allCond$cellType[int.allCond$cellType %in% c("3")] <- "HSC"
int.allCond$cellType[int.allCond$cellType %in% c("9")] <- "cholangiocytes"
int.allCond$cellType[int.allCond$cellType %in% c("10")] <- "T cells"
int.allCond$cellType[int.allCond$cellType %in% c("13")] <- "B cells"
int.allCond$cellType[int.allCond$cellType %in% c("15")] <- "DC/monocyte"
int.allCond$cellType[int.allCond$cellType %in% c("18")] <- "doublets"
int.allCond$cellType[int.allCond$cellType %in% c("17")] <- "proliferating cells"

int.allCond$cellType <- factor(int.allCond$cellType,levels = c("hepatocytes",
                                                               "KC",
                                                               "doublets",
                                                               "T cells",
                                                               "DC/monocyte",
                                                               "HSC",
                                                               "cholangiocytes",
                                                               "B cells",
                                                               "LSEC",
                                                               "proliferating cells"
))


int.allCond_woDoublets <- subset(x=int.allCond,subset = cellType %in% c("hepatocytes",
                                                                        "KC",
                                                                        #"doublets",
                                                                        "T cells",
                                                                        "DC/monocyte",
                                                                        "HSC",
                                                                        "cholangiocytes",
                                                                        "B cells",
                                                                        "LSEC",
                                                                        "proliferating cells"
))
Idents(int.allCond_woDoublets) <- "cellType"

int.allCond_woDoublets$cellType <- factor(int.allCond_woDoublets$cellType,levels = c("hepatocytes",
                                                                                     "KC",
                                                                                     "T cells",
                                                                                     
                                                                                     "HSC",
                                                                                     "cholangiocytes",
                                                                                     "B cells",
                                                                                     "LSEC",
                                                                                     "DC/monocyte",
                                                                                     "proliferating cells"
),ordered = T)

DimPlot(int.allCond_woDoublets, 
        reduction = "INTwnn.umap",
        label = TRUE,
        label.size = 5,
        repel = TRUE,
        split.by = "orig.ident") + 
  ggtitle("integrated WNN")

# Cleaning object to reduce environment size
saveRDS(int.allCond_woDoublets,"./int.allCond_woDoublets.RDS")
rm(list = setdiff(ls(), lsf.str()))

# Subclustering myeloid ----
int.allCond_woDoublets <- readRDS("./int.allCond_woDoublets.RDS")

int.allCond_myeloid <- subset(x=int.allCond_woDoublets,subset= cellType %in% c("KC"))
DefaultAssay(int.allCond_myeloid) <- "integrated"
Idents(int.allCond_myeloid) <- "cellType"
int.allCond_myeloid <- FindMultiModalNeighbors(
  int.allCond_myeloid, reduction.list = list("pca", "harmony"), 
  dims.list = list(1:50, 2:50),
  k.nn = 60)

int.allCond_myeloid <- RunUMAP(int.allCond_myeloid, 
                               nn.name = "weighted.nn", 
                               reduction.name = "umap.myeloid.wnn", 
                               reduction.key = "umap.myeloid.wnn", 
                               assay = "integrated")
int.allCond_myeloid <- FindClusters(int.allCond_myeloid, 
                                    graph.name = "wsnn", 
                                    algorithm = 3, 
                                    verbose = FALSE)


int.allCond_myeloid$merged <- paste0(int.allCond_myeloid$orig.ident,"_",int.allCond_myeloid$seurat_clusters)

int.allCond_myeloid$merged <- paste0("KC_",int.allCond_myeloid$seurat_clusters)


## cell number cutoff ----
# we take out cluster 5 as it has below <50 cells in total
table(int.allCond_myeloid$seurat_clusters)
int.allCond_myeloid_woCluster5 <- subset(int.allCond_myeloid, seurat_clusters!="5")

## After subclustering find markers ----

DefaultAssay(int.allCond_myeloid_woCluster5) <- "RNA"
int.allCond_myeloid_woCluster5 <- NormalizeData(int.allCond_myeloid_woCluster5)
MarkerList <- FindMarkers2Excel(
  int.allCond_myeloid_woCluster5,
  idents = "orig.ident",
  grouping.var = "seurat_clusters",
  filename = "./MarkerList_integratedWNN_myeloid.xlsx"
)

colorSubset<-c("#cb997e","#2a9d8f",
               "#8ab17d","#e9c46a",
               "#f4a261","#e76f51")
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

# CellChat analysis ----
colorTheme = c("#c6c6c6","#606060","#c12c38","#e0775f","#f3b694","#fce2d0")
names(colorTheme) <- c("wtCDCDCD",
                       "wtCDCDHFD",
                       "wtHFDCDCD",
                       "wtHFDCDHFD",
                       "wtHFDHFDCD",
                       "wtHFDHFDHFD")

DefaultAssay(int.allCond_woDoublets) <- "RNA"

# identify KC_5 cells to remove
removeCells <- names(int.allCond_myeloid$seurat_clusters)[int.allCond_myeloid$seurat_clusters == "5"] 
allToInclude <- setdiff(names(int.allCond_woDoublets$seurat_clusters),removeCells)
int.allCond_woDoublets_woclust5 <- int.allCond_woDoublets[,allToInclude]


int.allCond_woDoublets_wtCDCDCD_KCsub <- subset(int.allCond_woDoublets_woclust5,orig.ident == "wtCDCDCD")
int.allCond_woDoublets_wtCDCDCD_KCsub <- NormalizeData(int.allCond_woDoublets_wtCDCDCD_KCsub)
int.allCond_woDoublets_wtHFDCDCD_KCsub <- subset(int.allCond_woDoublets_woclust5,orig.ident == "wtHFDCDCD")  
int.allCond_woDoublets_wtHFDCDCD_KCsub <- NormalizeData(int.allCond_woDoublets_wtHFDCDCD_KCsub)

## wtCDCDCD prep ----
data.input <- int.allCond_woDoublets_wtCDCDCD_KCsub[["RNA"]]@data # normalized data matrix
int.allCond_woDoublets_wtCDCDCD_KCsub$samples <- int.allCond_woDoublets_wtCDCDCD_KCsub$orig.ident
table(Idents(int.allCond_woDoublets_wtCDCDCD_KCsub)) # ensure this are the cell labels
labels <- Idents(int.allCond_woDoublets_wtCDCDCD_KCsub) 


meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
meta$labels <- as.character(meta$labels)

Idents(int.allCond_myeloid) <- int.allCond_myeloid$merged
meta[names(Idents(int.allCond_myeloid)[int.allCond_myeloid$orig.ident == "wtCDCDCD"]),"labels"] <- paste0(Idents(int.allCond_myeloid)[int.allCond_myeloid$orig.ident == "wtCDCDCD"])

Idents(int.allCond_woDoublets_wtCDCDCD_KCsub) <- as.factor(meta$labels)
Idents(int.allCond_woDoublets_wtCDCDCD_KCsub) <- gsub("wtCDCDCD","KC",Idents(int.allCond_woDoublets_wtCDCDCD_KCsub))

cellChat_wtCDCDCD <- createCellChat(object = int.allCond_woDoublets_wtCDCDCD_KCsub, 
                                    #meta = meta,
                                    group.by = "ident", 
                                    assay = "RNA")

## wtHFDCDCD prep ----
data.input <- int.allCond_woDoublets_wtHFDCDCD_KCsub[["RNA"]]@data # normalized data matrix
int.allCond_woDoublets_wtHFDCDCD_KCsub$samples <- int.allCond_woDoublets_wtHFDCDCD_KCsub$orig.ident

labels <- Idents(int.allCond_woDoublets_wtHFDCDCD_KCsub)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
meta$labels <- as.character(meta$labels)
meta[names(Idents(int.allCond_myeloid)[int.allCond_myeloid$orig.ident=="wtHFDCDCD"]),"labels"] <- paste0(Idents(int.allCond_myeloid)[int.allCond_myeloid$orig.ident=="wtHFDCDCD"])

Idents(int.allCond_woDoublets_wtHFDCDCD_KCsub) <- as.factor(meta$labels)
Idents(int.allCond_woDoublets_wtHFDCDCD_KCsub) <- gsub("wtHFDCDCD","KC",Idents(int.allCond_woDoublets_wtHFDCDCD_KCsub))

cellChat_wtHFDCDCD <- createCellChat(object = int.allCond_woDoublets_wtHFDCDCD_KCsub, 
                                     group.by = "ident", 
                                     assay = "RNA")

## prep database ----
CellChatDB <- CellChatDB.mouse 
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB,non_protein = F)
cellChat_wtHFDCDCD@DB <- CellChatDB.use
cellChat_wtCDCDCD@DB <- CellChatDB.use

# Note that one can use projected data (projection of gene expression value to their neighbours within
# a PP-Network)

# subset the expression data of signaling genes for saving computation cost
cellChat_wtCDCDCD <- subsetData(cellChat_wtCDCDCD) # This step is necessary even if using the whole database
cellChat_wtHFDCDCD <- subsetData(cellChat_wtHFDCDCD)
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 2 * 1024^3)  # Increase to 2 GB

## CellChat ----
cellChat_wtCDCDCD <- identifyOverExpressedGenes(cellChat_wtCDCDCD)
cellChat_wtCDCDCD <- identifyOverExpressedInteractions(cellChat_wtCDCDCD)

cellChat_wtHFDCDCD <- identifyOverExpressedGenes(cellChat_wtHFDCDCD)
cellChat_wtHFDCDCD <- identifyOverExpressedInteractions(cellChat_wtHFDCDCD)

cellChat_wtCDCDCD <- computeCommunProb(cellChat_wtCDCDCD, type = "truncatedMean",trim = 0.1)
cellChat_wtHFDCDCD <- computeCommunProb(cellChat_wtHFDCDCD, type = "truncatedMean",trim = 0.1)

cellChat_wtCDCDCD <- filterCommunication(cellChat_wtCDCDCD, min.cells = 10)
cellChat_wtHFDCDCD <- filterCommunication(cellChat_wtHFDCDCD, min.cells = 10)

cellChat_wtCDCDCD <- computeCommunProbPathway(cellChat_wtCDCDCD)
cellChat_wtHFDCDCD <- computeCommunProbPathway(cellChat_wtHFDCDCD)

cellChat_wtCDCDCD <- aggregateNet(cellChat_wtCDCDCD)
cellChat_wtHFDCDCD <- aggregateNet(cellChat_wtHFDCDCD)


object.list_KCsub <- list(wtCDCDCD = cellChat_wtCDCDCD, wtHFDCDCD = cellChat_wtHFDCDCD)
saveRDS(object.list_KCsub, "object.list_KCsub_woCLust5.rds")

## Visualise results ----
cellchat_KCsub <- mergeCellChat(object.list_KCsub, add.names = names(object.list_KCsub))

gg1 <- rankNet(cellchat_KCsub, 
               mode = "comparison", 
               measure = "weight",
               sources.use = c("KC_0","KC_1","KC_2",
                               "KC_3","KC_4"), 
               targets.use = NULL, 
               stacked = T, do.stat = T,
               bar.w = 0.75,
               font.size = 12,
               #color.use = c("grey","red","green"),
               color.use = colorTheme[c(1,3)]
)

svglite::svglite(filename="rankNet_KCsub_toAll_woKC5.svg", 
                 width = 10, height = 20)
gg1
dev.off()

# TF analysis with decouplR ----
Idents(int.allCond_myeloid_woCluster5)<- int.allCond_myeloid_woCluster5$seurat_clusters
int.allCond_myeloid_woCluster5_SCTprepped <- PrepSCTFindMarkers(int.allCond_myeloid_woCluster5)
markers_loosened <- FindAllMarkers(int.allCond_myeloid_woCluster5_SCTprepped,
                                   min.pct = 0.10,
                                   #min.diff.pct=0.20,
                                   logfc.threshold = 0.20,
                                   return.thresh = 0.01,
                                   only.pos = F)
markersList <- list()

for(i in 0:4){
  markersList[[paste0("cluster",i)]] <- markers_loosened[markers_loosened$cluster==i,]
}
names(markersList)
lapply(markersList,nrow)

### Write to Excel

write.xlsx(markersList, file = "MarkerExpression_perKC_subcluster.xlsx")

allDE_pAdj <- lapply(markersList, function(x) x <- x[x$p_val_adj<0.05,])
allDE_pAdj_onlyGenes <- unique(unlist(lapply(allDE_pAdj, function(x) x$gene)))

mat <- as.matrix(int.allCond_myeloid_woCluster5@assays$SCT@data)

# cluster_wo5 <- int.allCond_myeloid_preppedSCT
cluster_wo5 <- int.allCond_myeloid_woCluster5
cluster_wo5$merged <- as.factor(paste0(cluster_wo5$orig.ident,"_",cluster_wo5$seurat_clusters))

mat_DEGenes <- mat[allDE_pAdj_onlyGenes,]
dim(mat_DEGenes)

net <- get_collectri(organism='10090', split_complexes=FALSE)

acts_DE_wo5 <- run_ulm(mat = mat_DEGenes, 
                       net=net, .source='source', .target='target',
                       .mor='mor', minsize = 10) 

cluster_wo5[['tfsulm_DE']] <- acts_DE_wo5 %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = cluster_wo5) <- "tfsulm_DE"

# Scale the data
cluster_wo5 <- ScaleData(cluster_wo5)
cluster_wo5@assays$tfsulm_DE@data <- cluster_wo5@assays$tfsulm_DE@scale.data


# Extract activities from object as a long dataframe

df <- t(as.matrix(cluster_wo5@assays$tfsulm_DE@data)) %>%
  as.data.frame() %>%
  mutate(cluster = cluster_wo5$merged) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score, trim = 0.25))

# Get top tfs with more variable means across clusters

sorted <- df %>%
  group_by(source) %>%
  summarise(std = max(mean)) %>%
  arrange(-abs(std))
sorted$Rank <- 1:nrow(sorted)
sorted$genes <- toupper(sorted$source)

n_tfs = nrow(sorted)
tfs <- df %>%
  group_by(source) %>%
  summarise(std = max(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

grepl("Hif1a",tfs)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

min(top_acts_mat)
max(top_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("darkblue", "white","red"))(palette_length)


my_breaks <- c(seq(-1, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 1, length.out=floor(palette_length/2)))


# annotate rows according to their group (HFD and CD)

anno_df <- data.frame(group=rownames(top_acts_mat))
rownames(anno_df) <- anno_df$group
anno_df$group <- gsub("_.*$","",anno_df$group)
anno_df$group  <- as.factor(anno_df$group)
ann_colors = list(group=c(wtCDCDCD="#c6c6c6",wtHFDCDCD="#c12c38"))

# Plot
library(pheatmap)


svglite::svglite("TF_analysis.svg")
pheatmap(top_acts_mat, 
         border_color = NA, 
         color=my_color,
         breaks = my_breaks,
         scale="none",
         cutree_rows=3,cutree_cols=3,
         cellwidth = 10, cellheight = 10,
         annotation_row = anno_df,
         annotation_colors = ann_colors
)
dev.off()


# Coverage plot ----
int.allCond_myeloid_woCluster5$merged_cond <- paste0(int.allCond_myeloid_woCluster5$orig.ident,
                                                     "_",
                                                     int.allCond_myeloid_woCluster5$seurat_clusters)
Idents(int.allCond_myeloid_woCluster5) <- int.allCond_myeloid_woCluster5$merged_cond
DefaultAssay(int.allCond_myeloid_woCluster5) <- "int_peaks"
# check tabixFiles
sample_typ <- "wtCDCDCD"
int.allCond_myeloid_woCluster5[["int_peaks"]]@fragments[[1]]@path <- paste0("/home/rstudio/program/data/",sample_typ,"/atac_fragments.tsv.gz")
sample_typ <- "wtHFDCDCD"
int.allCond_myeloid_woCluster5[["int_peaks"]]@fragments[[2]]@path<- paste0("/home/rstudio/program/data/",sample_typ,"/atac_fragments.tsv.gz")

regionQuery <- "chr7-19695000-19700000" # Apoe
#regionQuery <- "chr9-46227000-46230466" # Apoa1
regionQuery <- "chr7-19699000-19700000" # Apoe Zoom

extend.downstream = 0
extend.upstream = 0

query_parts <- strsplit(regionQuery, "-")[[1]]
chr <- query_parts[1]
start <- as.numeric(query_parts[2]) - extend.upstream
end <- as.numeric(query_parts[3]) + extend.downstream

# Create the updated regionQuery string
extendede_query <- paste(chr, start, end, sep = "-")

#get Diff peaks
DefaultAssay(int.allCond_myeloid_woCluster5) <- 'int_peaks'
Idents(int.allCond_myeloid_woCluster5) <- int.allCond_myeloid_woCluster5$merged_cond

# da_peaks_all <- list()
# for(i in 0:4){
#   DE <- FindMarkers(int.allCond_myeloid_woCluster5,
#                     ident.1 = paste0("wtHFDCDCD_",i),
#                     ident.2 = paste0("wtCDCDCD_",i),
#                     test.use = 'LR',
#                     min.pct = 0.0,
#                     logfc.threshold = 0,
#                     latent.vars = 'nCount_peaks',
#                     only.pos =F) 
#   da_peaks_all[[paste0("cluster",i)]] <- DE
# }
# 
# da_peaks_all_de <- lapply(da_peaks_all, function(x){
#   subset(x, p_val_adj < 0.1)
# })


overlap_peaks <- Signac::findOverlaps(int.allCond_myeloid_woCluster5, 
                                      StringToGRanges(extendede_query),
                                      type = "any",
                                      ignore.strand = T)

overlapRanges <- StringToGRanges("chr7-19699500-19700000") # for Apoe
#overlapRanges <- StringToGRanges("chr9-46228000-46228500") # for Apoa1

overlapRanges$color <- "#c12c38"

colorVector <- c(rep(c("#c6c6c6","#c12c38"),6))
names(colorVector) <- c("0_wtCDCDCD","0_wtHFDCDCD",
                        "1_wtCDCDCD","1_wtHFDCDCD",
                        "2_wtCDCDCD","2_wtHFDCDCD",
                        "3_wtCDCDCD","3_wtHFDCDCD",
                        "4_wtCDCDCD","4_wtHFDCDCD")

Idents(int.allCond_myeloid_woCluster5) <- int.allCond_myeloid_woCluster5$merged
covPlot <- Signac::CoveragePlot(
  object = int.allCond_myeloid_woCluster5,
  group.by = 'orig.ident',
  split.by = 'seurat_clusters',
  #region = "chr7-19700275-19701843",
  region = extendede_query,
  region.highlight = overlapRanges,
  annotation = F,
  #features = geneOfInterest,
  peaks = F,
  #extend.upstream = extend.downstream,#
  #extend.downstream = extend.downstream_plot#,
  #upstream=all_df[all_df$gene_name==geneOfInterest,"distance"],
  #downstream=all_df[all_df$gene_name==geneOfInterest,"distance"]
) + scale_fill_manual(values = colorVector)

gene_plot <- AnnotationPlot(
  object = int.allCond_myeloid_woCluster5,
  region = extendede_query
  #region = queryGene_Df$query_upstream
) + theme_classic(base_size = 20)

gene_plot$layers[[length(gene_plot$layers)]]$aes_params$size <- 5


# svglite::svglite("Apoe_region.svg", width = 10, height = 10)
CombineTracks(
  plotlist = list(covPlot, gene_plot),
  heights = c(10, 1),
  widths = c(10)
)
#dev.off()

## Indicate Ppar binding site
# get sequence of region of interest
sequence <- getSeq(BSgenome.Mmusculus.UCSC.mm10, 
                   names = unlist(strsplit(regionQuery,"-"))[1], 
                   start = start, 
                   end = end)


# Download PPAR motif

url <- "http://jaspar.elixir.no/temp/20241125144402_JASPAR2024_individual_matrices_1646816_pfm.zip"
destfile <- "PFM_matrices.zip"  # Destination file name
download.file(url, destfile, mode = "wb")

unzip_dir <- "PFM_matrices"  # Directory to extract files
unzip(destfile, exdir = unzip_dir)
extracted_files <- list.files(unzip_dir, full.names = TRUE)

url <- "http://jaspar.elixir.no/temp/20241125144549_JASPAR2024_individual_matrices_1646816_jaspar.zip"
destfile <- "JASPAR_matrices.zip"  # Destination file name
download.file(url, destfile, mode = "wb")

unzip_dir <- "JASPAR_matrices"  # Directory to extract files
unzip(destfile, exdir = unzip_dir)
extracted_files <- list.files(unzip_dir, full.names = TRUE)

TFs_binbdingSites <- overlapRanges

for(i in list.files(unzip_dir,full.names = T)){
  pfm_file <- i
  print(pfm_file)
  pfm <- TFBSTools::readJASPARMatrix(pfm_file, matrixClass = "PFM")
  pfm_single <- pfm[[1]]
  
  siteset <- searchSeq(toPWM(pfm_single), 
                       sequence, 
                       seqname="seq1", 
                       min.score="90%", 
                       strand="*")
  
  df_sites <- as(siteset, "data.frame")
  # filter for those which are before Apoe TSS
  if(nrow(df_sites)==0){
    print(paste0("skip:",pfm_file))
    next
  }
  
  TF_binding_site <- StringToGRanges(paste0("chr7-",
                                               start+df_sites$start,
                                               "-",
                                               start+df_sites$end))
  TF_binding_site$color <- "black"
  TF_binding_site$ID <-pfm_single@ID
  TF_binding_site$name <-pfm_single@name

  
  covPlot <- Signac::CoveragePlot(
    object = int.allCond_myeloid_woCluster5,
    group.by = 'orig.ident',
    split.by = 'seurat_clusters',
    region = extendede_query,
    region.highlight = c(overlapRanges,TF_binding_site),
    annotation = F,
    peaks = F,
  ) + scale_fill_manual(values = colorVector)

  gene_plot <- AnnotationPlot(
    object = int.allCond_myeloid_woCluster5,
    region = extendede_query
  ) + theme_classic(base_size = 20)
  
  final <- CombineTracks(
    plotlist = list(covPlot, gene_plot),
    heights = c(10, 1),
    widths = c(10)
  )
  
  svglite::svglite(paste0(pfm_single@name,"_Apoa_region.svg"), width = 15, height = 15)
  plot(final)
  dev.off()
  
  TFs_binbdingSites <- c(TFs_binbdingSites,TF_binding_site)
}

overlaps <- findOverlaps(TFs_binbdingSites[2:length(TFs_binbdingSites)])


# Do footprinting analysis ----
# library(JASPAR2022)
# pwm <- getMatrixSet(
#   x = JASPAR2022,
#   opts = list(species = 10090, all_versions = T)
# )
unzip_dir <- "JASPAR_matrices"
for(i in list.files(unzip_dir,full.names = T)){
  pfm_file <- i
  print(pfm_file)
  pfm <- TFBSTools::readJASPARMatrix(pfm_file, matrixClass = "PFM")
  TF_name <- names(pfm)

  
  # add motif information
  DefaultAssay(int.allCond_myeloid_woCluster5) <-"int_peaks"
  int.allCond_myeloid_woCluster5 <- AddMotifs(int.allCond_myeloid_woCluster5, 
                                              genome = BSgenome.Mmusculus.UCSC.mm10, 
                                              pfm = pfm)
  
  int.allCond_myeloid_woCluster5 <- Footprint(
    object = int.allCond_myeloid_woCluster5,
    motif.name = TF_name,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    in.peaks = F
  )
  
  
  
  
  Idents(int.allCond_myeloid_woCluster5) <- int.allCond_myeloid_woCluster5$orig.ident
  
  plottestFootprint <- PlotFootprint(int.allCond_myeloid_woCluster5, 
                features = TF_name)
  svglite::svglite(paste0(TF_name,"_tfFootprint.svg"))
  plot(plottestFootprint)
  dev.off()
  
}

  

PlotFootprint(int.allCond_myeloid_woCluster5, 
              features = c("Arnt"))
PlotFootprint(int.allCond_myeloid_woCluster5, 
              features = c("Pparg::Rxra"))

FoorPrintData <- GetFootprintData(int.allCond_myeloid_woCluster5,
                 features = "Arnt")

table(FoorPrintData$group)

subsetOfFootprint <- subset(FoorPrintData,class=="Observed")
# wtCDCDCD wtHFDCDCD 
# 506       506 
table(subsetOfFootprint$position)
head(subsetOfFootprint)

ggplot(subsetOfFootprint,aes(x=position,y=norm.value, color = group))+
  geom_point()

# subset to the different sets
wtCDCDCD_motif <- subset(int.allCond_myeloid_woCluster5,orig.ident=="wtCDCDCD")
wtHFDCDCD_motif <- subset(int.allCond_myeloid_woCluster5,orig.ident=="wtHFDCDCD")

peaks.with.motif_CD <- 
  GetMotifData(object = wtCDCDCD_motif)[, "MA0065.2"]
table(peaks.with.motif_CD)

peaks.with.motif_HFD <- 
  GetMotifData(object = wtHFDCDCD_motif)[, "MA0065.2"]
table(promoter.peaks.with.motif_HFD)
wtCDCDCD <- AccessiblePeaks(wtCDCDCD_motif)
wtHFDCDCD <- AccessiblePeaks(wtHFDCDCD_motif)

norm.data.CD <- GetAssayData(object = wtCDCDCD_motif, 
                          assay = 'peaks', 
                          slot = 'data')

norm.data.HFD <- GetAssayData(object = wtHFDCDCD_motif, 
                             assay = 'peaks', 
                             slot = 'data')

#try to plot heatmap
# y = position, y = binding site , fill Norm value, sep by orig ident.

motif_positions <- Motifs(int.allCond_myeloid_woCluster5)
arnt_motif <- "MA0497.1"  # Replace with the exact motif name or ID for Arnt

genome_ranges <- GRanges(seqnames = seqnames(BSgenome.Mmusculus.UCSC.mm10),
                         ranges = IRanges(start = 1, end = seqlengths(BSgenome.Mmusculus.UCSC.mm10)))

# Find matches for the Arnt motif across the genome
binding_sites <- motifmatchr::matchMotifs(
  pwms = pwm[[1]], 
  subject = genome_ranges, 
  genome = BSgenome.Mmusculus.UCSC.mm10, 
  out = "positions"  # Return positions of matches
)

# get all bindings sites
binding_sites_df <- as.data.frame(binding_sites)
head(binding_sites_df)
binding_sites_extended <- resize(binding_sites, width = 500, fix = "center")


# Extract peak values around binding sites for cells in conditions A and B
peak_matrix <- FeatureMatrix(
  fragments = Fragments(int.allCond_myeloid_woCluster5),
  features = binding_sites_extended[[1]],
  cells = Cells(int.allCond_myeloid_woCluster5)
)

# Split by condition
cells_wtCDCDCD <- Cells(int.allCond_myeloid_woCluster5)[int.allCond_myeloid_woCluster5$orig.ident == "wtCDCDCD"]
cells_wtHFDCDCD <- Cells(int.allCond_myeloid_woCluster5)[int.allCond_myeloid_woCluster5$orig.ident == "wtHFDCDCD"]

peaks_wtCDCDCD <- peak_matrix[, cells_wtCDCDCD]
peaks_wtHFDCDCD <- peak_matrix[, cells_wtHFDCDCD]


# Combine data into a format suitable for plotting
positions <- seq(-250, 250, length.out = nrow(binding_sites_df))
plot_data <- data.frame(
  Position = rep(positions, 2),
  Mean_Peak = c(rowSums(peaks_wtCDCDCD), rowSums(peaks_wtHFDCDCD)),
  Condition = rep(c("wtCDCDCD", "wtHFDCDCD)"), each = length(positions))
)

# Plot
hist(plot_data$Mean_Peak)
plot_data_sub <- plot_data[plot_data$Mean_Peak > 10,]

ggplot(plot_data_sub, aes(y = Mean_Peak, color = Condition)) +
  geom_boxplot() +
  #scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme_classic(base_size = 14) +
  labs(
    title = "TF Arnt Binding Sites Peak Intensities",
    x = "Position Relative to Motif (bp)",
    y = "Mean Peak Value",
    fill = "Condition"
  )

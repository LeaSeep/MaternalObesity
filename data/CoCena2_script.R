### R version check
if(grepl(c("3.4.4"),R.version.string)){
  print("Welcome to CoCena²")
}else{
  print("You can try but there can be troubles with other versions")
}
######################################################################

#####Set working directory
setwd("~/Documents/Paul_Creld2/Paul_v2")
mainDir <- getwd()
subDir <- c("CoCena_Drug_2")
dir.create(file.path(mainDir , subDir))
setwd(file.path(mainDir , subDir))
originalwd <- getwd()

# load packages

packages_to_load <- c("igraph" , "plotly" ,"ggplot2" ,"bnstruct", "gridExtra" , "Hmisc" , 
                      "RColorBrewer" , "gtools" , "rlist" , "reshape2" , "gplots" , "moduleColor" , 
                      "NMF" , "rlist" , "ggpubr","RJSONIO","httr","stringr","XML",
                      "devtools","clues","bnstruct","fields","curl","httpuv","intergraph","bnlearn","fields","viridis","foreign")

lapply(packages_to_load , require , character.only=TRUE)
install_github('cytoscape/cytoscape-automation/for-scripters/R/r2cytoscape')
install.packages("tictoc", dependencies=TRUE, type="source")
library(rJava)
library(r2cytoscape)
library(igraph)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(RColorBrewer)
library(gtools)
library(rlist)
library(reshape2)
library(gplots)
library(moduleColor)
library(NMF)
library(rlist)
library(ggpubr)
library(tictoc)
library(grid)
library(gridExtra)
library(viridis)
library(clues)
library(readr)
library(circlize)

library(foreign)
library(httpuv)

install.packages("jasmine")
library(devtools)
install_github("kassambara/r2excel")
#remove.packages(rJava)
library(r2excel)
library(grid)
library(clusterProfiler)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ReactomePA")
library(ReactomePA)
source("https://bioconductor.org/biocLite.R")
BiocManager::install("pcaGoPromoter")
library(pcaGoPromoter)
#source("https://bioconductor.org/biocLite.R")
BiocManager::install("biomaRt")
#biocLite("biomaRt")

library(biomaRt)
library(reshape2)
library(dplyr)

#####load in your data

###############################################################################################
#make sure that in your dataset, 
#gene names are as rownames and samples are as column names!
###############################################################################################
Dataset_1 <- batch_corrected_rld
rownames(Dataset_1)<-Dataset_1$ID
#rownames(Dataset_1)<-Dataset_1$ID
Dataset_1$ID<-NULL
#Dataset_1$TF_network_anno<-NULL
#Dataset_1$Tissue<-NULL
View(Dataset_1)
#Dataset_1<-t(Dataset_1)
###############################################################################################
#make sure that sample names are in rownames
#MAKE sure that there is no " - "in merged names - cytoscape cannot handle these
#Make sure that there is no space in your names- cytoscape cannot handle these
######################################################################################
info_Dataset <- annotation
rownames(info_Dataset)<-info_Dataset$ID
info_Dataset$ID<-NULL

View(info_Dataset)
info_Dataset$X.1 <-NULL

#####################################################################################
#set the organism your data is from
organism=c("mouse")
####################################################################################
##load additional needed Information!

library(clusterProfiler)
TF_list <- clusterProfiler::read.gmt("~/Downloads/c3.tft.v7.1.symbols.gmt.xls" , header = TRUE )
gmtfile <- clusterProfiler::read.gmt("~/Downloads/h.all.v6.2.symbols.gmt (1).txt")
epigenetic_modulators <- read_delim("Y:/Marie/oestreich.m/epigenetic_modulators.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
GO_immune_response <- xlsx.readFile("Y:/Marie/TheFinalProblem/GO_immune_response.xlsx", sheetIndex =1)
GO_metabolic_process <- xlsx.readFile("Y:/Marie/TheFinalProblem/GO_metabolic_process.xlsx", sheetIndex =1)
GO_signalling <- xlsx.readFile("Y:/Marie/TheFinalProblem/GO_signalling.xlsx", sheetIndex =1)
####################################################################################



#Save original data to new variable which could be changed during process
#original Dataset can always be found under Dataset_1
original_data <- Dataset_1


#####################################################################################
##filter data if wanted
filter=FALSE
#####################################################################################
#1) on transcription factors
original_data <- original_data[rownames(original_data) %in% TF_list$Human , ] ## HUman or mouse
#2) any other gene_list
gene_list<-UNION_DE
original_data <- original_data[rownames(original_data) %in% gene_list , ] ## HUman or mouse


####################################################################################
####################################################################################
collectGarbage()

#####
#first visualisation of your data
col.pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

original_data <- original_data[apply(original_data,1,function(x) !sd(x)==0),]

original_data <- original_data[1:500,]

pheatmap::pheatmap(original_data, 
                   cluster_row = F,
                   cluster_cols = F,
                   color = col.pal,
                   scale = c("row"),
                   annotation_col = info_Dataset,
                   main = c("Heatmap: Gene Expression Data - top 500 Dataset "))

####################################################################################
#There will be a function collecting all set variables values in this list "summary"
###
summary<-list()
####################################################################################
####################################################################################


###Start CoCena²
correlation_df<-correlation(type_of_correlation="pearson",            #pearson / spearman
                            pVal_treshold=0.05,                       #pValue Treshold for calculated correlation
                            save_raw_data=FALSE)                      #save correlation data in .txt

summary[["correlation"]]<-correlation_df$summary
correlation_df<-correlation_df$raw_data


####deciding on cutoff
library(igraph)
cutoff_options<-cutoff_visualisation(correlation_df=correlation_df,
                                     min_corr=0.75,                    # min correlation 
                                     range_cutoff_length=20,          #No.of cutoffs tested (all higher than min_cor)
                                     min_no_for_cluster=10)            #how many nodes are minimum for a network 

summary[["cutoff_visualisation"]]<-cutoff_options$summary

####################################################################################
##Look at plot
cutoff_options$Plotly_object
#slide through cutoffs and decide on one
#noch händisch einzutragen
####################################################################################
cutoff_options$Cutoff_df
chosen_cutoff<-0.791
###Change to you chosen cutoff!!
####################################################################################
cutoff_wd<-paste0(originalwd,"/",chosen_cutoff)
dir.create(file.path(cutoff_wd))
summary[["chosen_cutoff"]]<-chosen_cutoff
setwd(cutoff_wd)

##plot Network optional change layout
layout_options <- grep("^layout_" , ls("package:igraph") , value = TRUE)[-1]
layout_options
# [1] "layout_as_bipartite"  "layout_as_star"       "layout_as_tree"      
# [4] "layout_components"    "layout_in_circle"     "layout_nicely"       
# [7] "layout_on_grid"       "layout_on_sphere"     "layout_randomly"     
# [10] "layout_with_dh"       "layout_with_drl"      "layout_with_fr"      
# [13] "layout_with_gem"      "layout_with_graphopt" "layout_with_kk"      
# [16] "layout_with_lgl"      "layout_with_mds"      "layout_with_sugiyama"

#####################################################################################
#####test different layouts

layouts_on_list_coord<-test_layout(layouts_to_test=c("layout_with_fr",
                                                     "layout_with_kk",
                                                     "layout_with_lgl"),            #choose from layout_options
                                   min_nodes_number_for_network=10)

summary[["test_layout"]]<-layouts_on_list_coord[[c("summary")]]
layouts_on_list_coord<-layouts_on_list_coord$layouts_on_list_coord
dev.off()


####must run these command to get an igraph object

return_list<-plot_network(data=correlation_df,
                          layout=c("layout_with_lgl"),                          #choose anything       if you tested layouts change to layout_on_grid to fasten things up
                          min_nodes_number_for_network=20,
                          show_HC= TRUE)                                       #see HC after cutting whole Dataset
summary[["plot_network"]]<-return_list$summary
igraph_object<-return_list[["graph_object"]]

######Scale free Topology check
fit_power_law(igraph_object)

####layout can be chosen from
layout<-return_list[["layout"]]

####################################################################################
######if test_layout was run
names(layouts_on_list_coord)# names of your layout options
layout<-layouts_on_list_coord$layout_with_lgl
####################################################################################
#GFC calculation
GFC_all_genes<-GFC_calculation(normdata=original_data,
                               group=c("condition"),       #column name from info_Dataset
                               data_in_log=FALSE,          #is the inloaded data at some point put in log (e.g. DeSeq's rlog)
                               range_GFC=2.5)              #if calculated GFC value is over range_GFC it is 
#set to range_GFC value due to nice visulaisation
summary[["GFC_calculation"]]<-GFC_all_genes$summary
GFC_all_genes<-GFC_all_genes$GFC_all_genes
#####error !!!!!
plot_GFC_networks(igraph_object=igraph_object,
                  print_to_pdf=T,                      #GFC plots will be plotted in pdf - ATTENTION when further work will be done in Corel Draw not recommended
                  print_edges_png_nodes_pdf=c("none"))     #Options: each - every plot will be printed in one pdf(Nodes) and one png (edges)
#one - all plots will be in one pdf(Nodes) and one pmg(edges)
#none - nothing will be saved
dev.off()

##cluster algos
# cfg <- cluster_label_prop(g)
# cfg <- cluster_fast_greedy(g)
# cfg <- cluster_louvain(g)
# cfg<-cluster_infomap(g)
# cfg<-cluster_walktrap(g)
# cfg<-cluster_spinglass(g)
# cfg <- cluster_edge_betweenness(g) #ACHTUNG DAUERT EWIG
##

source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
library(ComplexHeatmap)
cluster_informationen<-heatmap_clustered(igraph_object=igraph_object,
                                         cluster_algo=c("cluster_louvain"),                                        #if set to auto optimal cluster algo will be chosen
                                         layout_for_network=layout,
                                         iterations=TRUE,                                               #iterations will produce a stable result but will take a little longer
                                         no_of_iterations=1,                                          #number of iterations
                                         max_cluster_count_per_gene=8,                                  #no of clusters which a gene is allowd to be advised to before putting it into waste cluster
                                         min_cluster_size=50,                                            # min size of cluster to be shown in heatmap
                                         desicion_to_see_plot=TRUE,                                     #whether to see the network plotted as well
                                         desicion_to_save_plot=FALSE,                                   #whether or not to save network 
                                         name_for_pdf_plot=c("Network_greater_clusters_clustered"),
                                         print_to_pdf=FALSE,                                           #whether or not to save heatmap
                                         name_for_pdf=c("Heatmap_greater_clusters_clustered"),
                                         average = "mean")

summary[["heatmap_clustered"]]<-cluster_informationen$summary
summary[["heatmap_clustered"]]
clustered_heatmap_data<-cluster_informationen$heatmap_df
clustered_heatmap_data
cluster_information<-cluster_informationen$color_cluster_min_size
#write.csv(cluster_informationen, "info.csv")
########
##unzufrieden mit der Heatmap?
##sie haben zu viele cluster die den selben trend über alle conditionen aufzeigen?
##sie möchten selber bestimmen wie viele cluster sie haben wollen ?
##dann benutzen sie folgende funktion
merged_cluster_data<-merge_cluster(clustered_heatmap_data,cluster_wanted = 5)
summary[["merge_cluster"]]<-merged_cluster_data$summary
cluster_information<-merged_cluster_data$cluster_information_new

library(pheatmap)
pheatmap(clustered_heatmap_data)


###################################################################################
#Circos Plot for all clusters

y_compare <- compareClusterGO() # This takes very long. This step can be skipped,
# then you just need to set the y_comp parameter
# in clusterCircos() to F.

heatmap_vec <- c("GFC_WT_Sucrose", "GFC_WT_Sucrose", "GFC_KO_TG", "GFC_WT_TG","GFC_WT_TM, GFC_KO_TM")

dev.off()

clusterCircos(hm_vec = TRUE, seg = 100, range = c(-1, 0, 1), y_comp = F)




################################to cytoscape########################################
###CYTOSCAPE MUST BE OPEN
toCytoscape_all(cluster_information=cluster_information)
####################################################################################

#testing out different cutoffs
search_for_good_cutoff(data=correlation_df,
                       min_nodes_number_for_network=10,            #min no of nodes connected to count as network
                       show_network=TRUE,                          #show plotted network
                       layout=layout_with_kk,                      #choose layout algorithm
                       cluster_algo= cluster_louvain,              # choose cluster algorithm
                       max_cluster_count_per_gene=8,               #no of clusters which a gene is allowd to be advised to before putting it into waste cluster
                       min_cluster_size=10,                         #min size of cluster to be shown in heatmap
                       abberation=0.1,                             #abberation from before chosen cutoff in percent , can be set to 0
                       no_cutoffs_tested=4,                        # no of cutoffs you want to test
                       cutoffs_to_test=c("0.5","0.8","0.2"))      #set abberation to 0 if you want certain cutoffs tested
####output ?! maybe list with the igraph objects ?!


unique(cluster_information$try)


############################################################################################################
##Cluster Profiler
##########################################################################################################
mart<-useMart("ensembl")
mart<-listDatasets(mart) # in dataframe names of available datasets attention only functions with human / mouse
human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") # load both independently from your organism
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl") # load both independently from your organism

clusterprofiler_results<-clusterprofiler_autoCena(cluster_to_check=c("all"),    #cluster to check can be either a certain cluster or all clusters
                                                  group=c("merged"))         #condition you wanna check must be a column name of your info_Dataset

summary[["clusterprofiler_results"]]<-clusterprofiler_results$summary
clusterprofiler_results$summary<-NULL
###
#following creates Dataframe where each gene gets its annotation and a color depending in what analysis it was called
# color_cluster_info<-cluster_information_to_color(vertex_attributes=cluster_information,        # Dataframe with every gene and its belonging cluster
#                                                  clusterprofiler_results=clusterprofiler_results) # cluster profiler results to add
# summary[["cluster_information_to_color"]]<-color_cluster_info$summary
# color_cluster_info<-color_cluster_info$color_saved
# color_cluster_info<-as.data.frame(color_cluster_info)

igraph_object_all<-igraph_object   # big network is saved in  igraph_object_all // parent network

#############################################################################################################
#############################################################################################################
########################################## Intraclusteral Analysis ##########################################
############################################################################################################
#cluster names you can chose from:
unique(cluster_information$try)

list_of_results<-plot_single_cluster(igraph_object=igraph_object_all,               #parent network 
                                     cluster_name=c("blue"),        #cluster you chosen
                                     top_percentage_for_hubs=0.25,                  #percentage of how many genes is allowed to be hubs
                                     allowed_edges=20,                              #allowed edges if gene is counted as hub
                                     allowed_edges_between_hubs=2,                  #allowed edges between hubs
                                     string_needs_to_be_redone=TRUE,
                                     string_treshold=500,                           #if testing different parameters you not always need to redo string -
                                     color_STRING=c("#c95555"),                     #color to mark edges stated in STRING
                                     color_edges=c("grey"),                         # color proposed edges /edges not found on STRING
                                     no_strings_to_be_string_hub=1,                   #String edges can bring forward new HUBS, 
                                     #determine how many you want to have allowed - 
                                     #percentage is depending on degree number
                                     label_all_TF=TRUE,                            #TRUE/FALSE is you want all Transcription factors labeld - 
                                     #if FALSE only labelled if these fall under categorie labelled
                                     color_label_if_TF=c("#D3436EFF"),              # color of TF
                                     color_label_normal=c("#231151FF"),             # color of not TF
                                     width_label_edge=2,                            #width label edge to be more striking
                                     width_normal=1,                                # width normal edge
                                     percentage_named=0.7,                          # how many genes you want to be labeld depending an all paramteres before
                                     size_label_boxes=30)                           # adjust label boxes size normally between 20 and 30
summary[["I-GIN"]]<-list_of_results$summary
to_save<-list_of_results$summary[2:length(list_of_results$summary)]
to_save<-list.rbind(to_save)
summary[["Output_IGIN"]]<-list_of_results[["igraph_object_small"]]
igraph_object_small<-list_of_results[["igraph_object_small"]]
vertex_attributes_label<-list_of_results[["vertex_attributes_cluster"]]
edges_selected_label<-list_of_results[["edges_attributes"]]
############################################################################################################
#########################################Plot this beautiful thing##########################################
###########################################################################################################
library(qgraph)

e <- igraph::get.edgelist(igraph_object_small,names = FALSE)

l <- qgraph::qgraph.layout.fruchtermanreingold(e,vcount=igraph::vcount(igraph_object_small),
                                               weights = igraph::edge.attributes(igraph_object_small)$weight,
                                               area=10*(igraph::vcount(igraph_object_small)^2*2),                          #play around with area and reoulse.rad do minimise overlap and
                                               repulse.rad=(igraph::vcount(igraph_object_small)^3.1),                        #achieve nice plotting
                                               niter = 1000)                                                             #the higher niter the longer it takes but also more trys to improve

#################
#IGIN
################
#pdf("Übelster_nicest_network.pdf", height = 10, width = 16)
#par(mar=c(5.1,4.1,4.1,2.1))

plot(igraph_object_small,layout=l,
     vertex.size=vertex_attributes_label$size,
     vertex.size2=vertex_attributes_label$size2,
     vertex.label=vertex_attributes_label$label,
     vertex.color=vertex_attributes_label$color_from_hierachy,
     vertex.label.color=vertex_attributes_label$label_color,
     vertex.frame.color=vertex_attributes_label$frame_color,
     edge.color=edges_selected_label$col,edge.curved=edges_selected_label$curved,
     vertex.label.font=2,vertex.shape=vertex_attributes_label$shape,
     edge.width=edges_selected_label$edge_width,
     main=to_save["cluster_name",])


#############################search and label certain genes by hand########################################
add_gene_label(gene_to_add=c("Trp53"))



##############################################################################
###plot bayesian_based intresting genes with GFC indicating it's expression###
bayesian_based_subnetworks(igraph_object_small)

pieplot_coloring_plot(igraph_object=igraph_object_small, 
                      degree_to_color=FALSE,degree_to_color_percentage=0.1,
                      bayesian_to_color=TRUE,
                      resize_nodes=TRUE,
                      size_parent_pies=15,
                      size_children_pies=10,
                      size_others=1,
                      edge_width=2)

dev.off()
# bayesian_network_igraph<-graph_from_data_frame(bnnet2$arcs,directed = TRUE)
# bayesian_network_igraph
# V(bayesian_network_igraph)$name
# label_options_bayesian<-data.frame(Gene=V(bayesian_network_igraph)$name, label = " ", color="lightgreen",size=4,stringsAsFactors = FALSE)
# str(label_options_bayesian)
# label_options_bayesian[label_options_bayesian$Gene %in% c("Mmp14"),"label"]<-c("Mmp14")
# label_options_bayesian[label_options_bayesian$Gene %in% c("Runx3"),"label"]<-c("Runx3")
# label_options_bayesian[label_options_bayesian$Gene %in% c("Nek6"),"label"]<-c("Nek6")
# label_options_bayesian[label_options_bayesian$Gene %in% c("Mmp14"),"color"]<-"yellow"
# label_options_bayesian[label_options_bayesian$Gene %in% c("Runx3"),"color"]<-"yellow"
# label_options_bayesian[label_options_bayesian$Gene %in% c("Nek6"),"color"]<-"yellow"
# label_options_bayesian[label_options_bayesian$Gene %in% c("Mmp14"),"size"]<-15
# label_options_bayesian[label_options_bayesian$Gene %in% c("Runx3"),"size"]<-15
# label_options_bayesian[label_options_bayesian$Gene %in% c("Nek6"),"size"]<-15
# shortestPath<-shortest_paths(bayesian_network_igraph,from="Arhgap24",to="Runx3")
# 
# E(bayesian_network_igraph, path = shortestPath$vpath[[1]])$color<-"red"
# graph_from_edgelist()
# 
# plot(bayesian_network_igraph,
#      layout=layout_as_tree,
#      #vertex.color=label_options_bayesian$color,
#      # vertex.label=label_options_bayesian$label,
#      vertex.size=5,
#      edge.arrow.size=0.5)


######saving
mainDir <- cutoff_wd
subDir <- c("IGIN")
dir.create(file.path(mainDir , subDir))
setwd(file.path(mainDir , subDir))
IGIN_wd<-getwd()



write.table(to_save,paste0(IGIN_wd,"/",to_save["cluster_name",],".txt"),sep="\t")



#####LEGENDE eher semi

# #par(mar=c(5.1,0.5,4.1,0.5))
# colfunc <- colorRampPalette(c('#B1B1B2','#982D80FF'))
# graphics::legend(x=-2,y=1.5, legend = c("Edges","known interactions","unknown interactions"," ","Labels","transcription factor","no transcription factor"," ","Nodes"),
#        col = c("white","#c95555","grey","white","white","#D3436EFF","#231151FF", "white") , pch = c(15,15,15,15,15,15,15,15),
#        bty = "n", pt.cex = 2.2, cex = c(1,1,1,1,1,1,1), horiz = F , text.font = c(2,1,1,1,2,1,1,2))
# 
# 
# xl <- 1
# yb <- 0.4
# xr <- 1.5
# yt <- -0.1
# 
# graphics::rect(-1.93,-head(seq(yb,yt,(yt-yb)/20),-1),-1.88,-tail(seq(yb,yt,(yt-yb)/20),-1), col = colfunc(20), border = NA)
# 
# 
# 
# 
# mtext(c("receiving > regulating"," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","regulating > receiving"),side=4,at=-tail(seq(yb,yt,(yt-yb)/10),20),line=-62,las=2,cex=1)
# 
# 
# 
# dev.off()



#############################################################################################################
#################################################  summary  ##################################################
##############################################################################################################
#big list with included data
summary #can load some time !

#output most important set values during script
summarized_summary<-summary_wrapped_up(print_text = TRUE)

list.save(summary, "summary.rdata") # will save all data which is important to run the script, there is another script
# optimised for this summary data which can be passed along and reloaded on every other computer
#list.load("summary.rdata")


############################################################################################################
############################################################################################################
########################################### MARIE WORK #####################################################
############################################################################################################

BayesOutput <- BayesNet(clust.name =c("turquoise") ,
                        cutoff = chosen_cutoff , 
                        CPratio = 5 ,
                        PCratio = 3)

BayesOutput$BayesCompleteList$Mmp14.children
BayesOutput$BayesCompleteList$Mmp14.parents
BayesOutput$BayesCompleteList$Runx3.children
BayesOutput$BayesCompleteList$Runx3.parents
#################################################################
##testing genes

# test_genes<-c("Runx3")
# test_genes_expression<-Dataset_1[test_genes,]
# test_genes_expression<-rbind(test_genes_expression,as.character(info_Dataset$merged))
# test_genes_expression<-rbind(test_genes_expression,ifelse(grepl("TAM", info_Dataset$merged),c("TAM"),c("Control")))
# test_genes_expression<-as.data.frame(t(test_genes_expression))
# test_genes_expression$gene<-test_genes
# colnames(test_genes_expression)<-c("expression","Condition","group","Gene")
# # condition_tam<-data.frame(test_genes_expression[,grepl("TAM",info_Dataset$merged)],condition="condition",check.names=FALSE)
# #
# # control<-data.frame(test_genes_expression[,!(grepl("TAM",info_Dataset$merged))],condition="control",check.names=FALSE)
# test_genes_expression$expression<-as.numeric(as.character(test_genes_expression$expression))
# 
# test_genes_expression$group<-as.factor(test_genes_expression$group)
# test_genes_expression$Condition<-as.factor(test_genes_expression$Condition)
# Mmp14_test_genes_expression<-test_genes_expression
# test_genes<-c("Mmp14")
# test_genes_expression<-Dataset_1[test_genes,]
# test_genes_expression<-rbind(test_genes_expression,as.character(info_Dataset$merged))
# test_genes_expression<-rbind(test_genes_expression,ifelse(grepl("TAM", info_Dataset$merged),c("TAM"),c("Control")))
# test_genes_expression<-as.data.frame(t(test_genes_expression))
# test_genes_expression$gene<-test_genes
# colnames(test_genes_expression)<-c("expression","Condition","group","Gene")
# # condition_tam<-data.frame(test_genes_expression[,grepl("TAM",info_Dataset$merged)],condition="condition",check.names=FALSE)
# #
# # control<-data.frame(test_genes_expression[,!(grepl("TAM",info_Dataset$merged))],condition="control",check.names=FALSE)
# test_genes_expression$expression<-as.numeric(as.character(test_genes_expression$expression))
# 
# test_genes_expression$group<-as.factor(test_genes_expression$group)
# test_genes_expression$Condition<-as.factor(test_genes_expression$Condition)
# Runx3_test_genes_expression<-test_genes_expression
# 
# all<-rbind(Runx3_test_genes_expression,Mmp14_test_genes_expression)
# all$Gene<-as.factor(all$Gene)
# 
# all_without_PRE<-all[!(grepl("PRE",all$Condition)),]
# 
# ggplot(all,aes(x=group,y=expression,fill=Gene))+geom_boxplot()+geom_jitter()
# 
# dev.off()

#######
Dataset1
blue <- cluster_information$Gene[cluster_information$try=="blue"]
blue_counts <- original_data[tolower(rownames(original_data)) %in% tolower(blue), ]
write.csv(skyblue2_counts, file="skyblue2_counts_all.csv")

darkorange <- cluster_information$Gene[cluster_information$try=="darkorange"]
darkorange_counts <- original_data[tolower(rownames(original_data)) %in% tolower(darkorange), ]
write.csv(darkorange_counts, file="darkorange_counts_all.csv")

antiquewhite4 <- cluster_information$Gene[cluster_information$try=="antiquewhite4"]
antiquewhite4_counts <- original_data[tolower(rownames(original_data)) %in% tolower(antiquewhite4), ]
write.csv(antiquewhite4_counts, file="antiquewhite4_counts_all.csv")

brown4 <- cluster_information$Gene[cluster_information$try=="brown4"]
brown4_counts <- original_data[tolower(rownames(original_data)) %in% tolower(brown4), ]
write.csv(brown4_counts, file="brown4_counts_all.csv")

mediumorchid <- cluster_information$Gene[cluster_information$try=="mediumorchid"]
mediumorchid_counts <- original_data[tolower(rownames(original_data)) %in% tolower(mediumorchid), ]
write.csv(mediumorchid_counts, file="mediumorchid_counts_all.csv")

royalblue <- cluster_information$Gene[cluster_information$try=="royalblue"]
royalblue_counts <- original_data[tolower(rownames(original_data)) %in% tolower(royalblue), ]
write.csv(royalblue_counts, file="royalblue_counts_all.csv")

honeydew1 <- cluster_information$Gene[cluster_information$try=="honeydew1"]
honeydew1_counts <- original_data[tolower(rownames(original_data)) %in% tolower(honeydew1), ]
write.csv(honeydew1_counts, file="honeydew1_counts_all.csv")



library(DOSE)

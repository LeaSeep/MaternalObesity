## create expert knowledge dot-plot based on ORA results
# Requires an xlsx sheet of the selected terms (can be found in data folder)
setwd("transcriptome_analysis")
library(clusterProfiler)
ORA_all_results <- readRDS("Transcriptomics_ORA_CoCena_results_KC_WT.rds")
terms_selected <- as.data.frame(readxl::read_excel("../data/terms for figure 2.xlsx",sheet = "pickedTermsOfInterest"))


# data selection
selectedTerms <- list()
for(i in names(ORA_all_results)){
  # intermeditate!!!
  selectedTerms_cluster <- subset(terms_selected,ClusterName==i)
  data_cluster = ORA_all_results[[i]]
  data_term=c()
  
  if(any(grepl(c("KEGG"),selectedTerms_cluster$Term))){
    data_tmp=ORA_all_results[[i]]$KEGG

    data_term <- rbind(data_term,
                   data_tmp[data_tmp$Description%in%toupper(selectedTerms_cluster$Term),])
  }
  if(any(grepl(c("GO"),selectedTerms_cluster$Term))){
    terms_selected_tmp <- trimws(gsub("GO_","",selectedTerms_cluster$Term))
    data_tmp=ORA_all_results[[i]]$GO
    data_term <- rbind(data_term,
                   data_tmp[data_tmp$Description%in%terms_selected_tmp,])
  }
  if(any(grepl(c("HALLMARK"),selectedTerms_cluster$Term))){
    data_tmp=ORA_all_results[[i]]$HALLMARK
    data_term <- rbind(data_term,
                   data_tmp[data_tmp$Description%in%selectedTerms_cluster$Term,])
  }
  data_term$type <- i
  selectedTerms[[i]] <- data_term

}

data <- rlist::list.rbind(selectedTerms)


library(ggplot2)

data$roundedGeneRatio <- 0

for(i in 1:nrow(data)){
  tmp <- as.numeric(unlist(strsplit(data[i,"GeneRatio"],split ="\\/")))
  data[i,"roundedGeneRatio"] <- round(tmp[1]/tmp[2],2)
}

data$type<-factor(data$type, levels = c("plum","steelblue","gold","turquoise","lightgreen" ),ordered = T)
data$Description <- gsub("_"," ",data$Description)

ggplot_v1<-
  ggplot(data)+
  geom_point(aes(x=roundedGeneRatio, y=Description,color = p.adjust),size=4)+
  theme_bw()+
  facet_grid(type~.,switch ="y",scales = "free_y")+
  theme(strip.placement = "outside")+
  ylab("")+
  xlab("Gene Ratio")+ 
  scale_x_continuous(breaks=c(0,0.1,0.2,0.3))+
  scale_color_gradient(low = "red", high = "blue")+ 
  theme(aspect.ratio=1)


g <- ggplot_gtable(ggplot_build(ggplot_v1))
stripr <- which(grepl('strip-l', g$layout$name))
fills <- levels(data$type)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)

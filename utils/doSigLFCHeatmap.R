# Parse your data like the following sample table
# data:
# entities x samples

# annotationEntities:
# entities x groups (e.g. Class)
# summarise_by = character (must be one of the colnames of annotationEntities)
# sampleAnno:

# samples x groups (e.g. wt or ko)
# LFC_between = character (must be one of the colnames of sampleAnno)
# FC_ctrl = character must be a level of LFC_between specified column
# givenFilename = for saving the heatmap (incl. extension to specify filetype)

# Outputs: 
#   - Heatmap
#   - stat analysis results


doSigLFCHeatmap=function(
    data_matrix,
    annotation,
    sample_anno_lipid,
    summarise_by="CLASS",
    LFC_between="diet",
    FC_ctrl="CDCDCD",
    givenFilename="Heatmap.png",
    subset = NULL,
    colorTheme,
    transposed = T
){

  paste0("Number of Lipids with negative entries (removed): ",length(unique(as.data.frame(which(data_matrix<0,arr.ind = T))[,1])))
  if(length(unique(as.data.frame(which(data_matrix<0,arr.ind = T))[,1]))>0){
    df_tmp=data_matrix[-which((data_matrix < 0), arr.ind = T)[,1],]
  }else{
    df_tmp=data_matrix
  }

  # Make mean per summarise_by group for each sample ----
  df_tmp=as.data.frame(cbind(df_tmp,annotation[rownames(df_tmp),summarise_by]))
  colnames(df_tmp)[ncol(df_tmp)] = summarise_by
  df_tmp[,-ncol(df_tmp)] = apply(df_tmp[,-ncol(df_tmp)],2,as.numeric)
  mean_per_lipid_per_condition = aggregate(x = df_tmp[,-ncol(df_tmp)],
                                           by = list(as.factor(df_tmp[,summarise_by])),
                                           FUN = mean)
  rownames(mean_per_lipid_per_condition) = mean_per_lipid_per_condition$Group.1
  mean_per_lipid_per_condition$Group.1 = NULL
  mean_per_lipid_per_condition_t = as.data.frame(t(mean_per_lipid_per_condition))
  df_tmp=cbind(mean_per_lipid_per_condition_t,sample_anno_lipid[rownames(mean_per_lipid_per_condition_t),])
  
  mean_per_lipid_class = aggregate(x = df_tmp[,1:length(unique(annotation[,summarise_by]))],     
                                   
                                   # Specify group indicator
                                   by = list(as.factor(df_tmp[,LFC_between])),      
                                   
                                   # Specify function (i.e. mean)
                                   FUN = mean)
  
  rownames(mean_per_lipid_class) = mean_per_lipid_class$Group.1
  mean_per_lipid_class$Group.1 = NULL
  
  # get FC's
  FC_table_per_lipidClass=apply(mean_per_lipid_class,2,function(x){
    x[setdiff(unique(df_tmp[,LFC_between]),FC_ctrl)]/x[FC_ctrl]
  })
  #get LFC
  LFC_table_per_lipidClass=log2(FC_table_per_lipidClass)
  
  #get sig
  result_sig = LFC_table_per_lipidClass
  result_sig[1:nrow(LFC_table_per_lipidClass),1:ncol(LFC_table_per_lipidClass)]=""
  
  # get Data frame of real p-vals
  result_sig_pVals = LFC_table_per_lipidClass
  result_sig_pVals[1:nrow(LFC_table_per_lipidClass),1:ncol(LFC_table_per_lipidClass)]=NA
    
  # prep factors
  df_tmp[,LFC_between] = as.factor(df_tmp[,LFC_between])
  df_tmp[,LFC_between] <- relevel(df_tmp[,LFC_between], ref = FC_ctrl)
  
  for(i in colnames(result_sig)){
    formula_to_test = paste0("`",i,"`","~ ",LFC_between)
    fit <- lm(formula = formula_to_test, data = df_tmp)
    result_testing=summary(fit) # to get the mean now intercept + estimate
    rownames(result_testing$coefficients)=gsub(LFC_between,"",rownames(result_testing$coefficients))
    result_sig_pVals[names(result_testing$coefficients[-1,4]),i] <- result_testing$coefficients[-1,4]
    is_one_star=result_testing$coefficients[-1,4]<0.05 # without intercept
    is_two_star=result_testing$coefficients[-1,4]<0.01 
    is_three_star=result_testing$coefficients[-1,4]<0.001
    result_sig[names(is_one_star),i][is_one_star]="*"
    result_sig[names(is_two_star),i][is_two_star]="**"
    result_sig[names(is_three_star),i][is_three_star]="***"
  }
  #plotting
  paletteLength=19 # needs to be uneven 

  # subset if wanted
  if(!is.null(subset)){
    LFC_table_per_lipidClass=LFC_table_per_lipidClass[,subset]
    result_sig=result_sig[,subset]
    result_sig_pVals=result_sig_pVals[,subset]
  }
  
  # Sometime small LFC but highly significant => adjust color
  onlySigLFC <- LFC_table_per_lipidClass[result_sig_pVals<0.05]
  
  if(any(onlySigLFC > 0)==T){
    min_positive_sig_LFC <- min(onlySigLFC[onlySigLFC>0])
  }else{
    min_positive_sig_LFC <- quantile(LFC_table_per_lipidClass[LFC_table_per_lipidClass>0],probs = 0.25)
  }

  if(any(onlySigLFC <0 )==T){
    min_negative_sig_LFC <- max(onlySigLFC[onlySigLFC<0])
  }else{
    min_negative_sig_LFC  <- quantile(LFC_table_per_lipidClass[LFC_table_per_lipidClass<0],probs = 0.75)
  }

  min_both <- min(abs(c(min_positive_sig_LFC,min_negative_sig_LFC)))

  
  
  nBins_each = floor(max(LFC_table_per_lipidClass)/min_both)
  
  myColor <- c(colorRampPalette(c("#0077b6", "#d9effa"))(nBins_each),
               "white",
               colorRampPalette(c("#ffe8e8","#D62828"))(nBins_each)
  )
  
  
  # if scaling wanted
  #LFC_table_per_lipidClass=scale(LFC_table_per_lipidClass)

  nbins_needed = ceiling((max(LFC_table_per_lipidClass)-min_both)/min_both)*2+1
  
  myBreaks <- (c(rev(seq(1,nbins_needed,by=2))*-1,seq(1,nbins_needed,by=2))*min_both)/2
  
  
  if(length(myBreaks)>20){
    myColor <- c(colorRampPalette(c("#0077b6", "#d9effa"))(10),
                 "white",
                 colorRampPalette(c("#ffe8e8","#D62828"))(10)
    )
    myBreaks <- c(
      seq(-max(LFC_table_per_lipidClass),-min_both,length.out=10),
      (min_both/-2),
      min_both/2,
      seq(min_both,max(LFC_table_per_lipidClass),length.out=10)
    )
  }
 
  ## Add group colors
  # browser()
  annoCol <- list(group = colorTheme)
  col_anno = data.frame(group = names(colorTheme))
  rownames(col_anno) = col_anno$group
  LFC_table_per_lipidClass = as.data.frame(LFC_table_per_lipidClass)
  if(transposed){
    data2PLot <- as.data.frame(t(LFC_table_per_lipidClass))
    resultSig2Plot <- as.data.frame(t(result_sig))
    col_2_row <- NULL
    col_2_column <- col_anno
    cluster_rows <- T
    cluster_cols <- T
  }else{
    data2PLot <- as.data.frame((LFC_table_per_lipidClass))
    resultSig2Plot <- as.data.frame((result_sig))
    col_2_row <- col_anno
    col_2_column <- NULL
    cluster_rows <- T
    cluster_cols <- T
  }
  
  heatmap <- pheatmap(data2PLot,
                   color = myColor,
                   breaks = myBreaks, # no manual adjustment needed if scaling 
                   display_numbers = resultSig2Plot, # needs to be a matrix of same dim as LFC Table
                   fontsize_number = 15,
                   cellheight = 15,
                   cellwidth = 20,
                   filename = givenFilename,
                   annotation_colors = annoCol,
                   annotation_row = col_2_row,
                   annotation_col = col_2_column,
                   cluster_rows = cluster_rows,
                   cluster_cols = cluster_cols,
  )
  # save as svg
  heatmap <- pheatmap(data2PLot,
                      color = myColor,
                      breaks = myBreaks, # no manual adjustment needed if scaling 
                      display_numbers = resultSig2Plot, # needs to be a matrix of same dim as LFC Table
                      fontsize_number = 15,
                      cellheight = 15,
                      cellwidth = 20,
                      filename = gsub(".png",".tiff",givenFilename),
                      annotation_colors = annoCol,
                      annotation_row = col_2_row,
                      annotation_col = col_2_column,
                      cluster_rows = cluster_rows,
                      cluster_cols = cluster_cols,
  )
  

  
  # create color histogram
  col_hist_df=data.frame(dataPoints=unlist(LFC_table_per_lipidClass))
  col_hist_df$breaks = cut(col_hist_df$dataPoints,breaks = myBreaks)
  
  color_df = data.frame(x=levels(col_hist_df$breaks),y=1)
  #color_df = rbind(df_for_color,df_for_color)

  
  library(ggplot2)
  color_df$x=factor(color_df$x,levels= levels(col_hist_df$breaks),ordered = T)
  #col_hist_df$breaks=factor(col_hist_df$breaks,levels= levels(color_df$breaks),ordered = T)
  #max(table(col_hist_df$breaks))
  col_hist_df$breaks <- as.factor(col_hist_df$breaks)
  

  colorKey=ggplot()+
    geom_raster(data=color_df,
                aes(x = x, y = 1,fill = x),
                alpha=0.8, show.legend = FALSE,
                vjust = 0)+
    geom_bar(data = col_hist_df,
             aes(x=breaks, y=(..count..)/(max(table(col_hist_df$breaks))+1)),
                   stat="count",
                   fill=NA,
                   color="black",
                   show.legend = FALSE)+
    scale_fill_manual(values = myColor,breaks = levels(color_df$x), drop=F)+
    scale_x_discrete(drop=FALSE,limits = levels(color_df$x))+
    theme_void()
  
  
  
  ggsave(filename = gsub(".png","_colorHist.png",givenFilename), plot = colorKey)
  ggsave(filename = gsub(".png","_colorHist.svg",givenFilename), plot = colorKey)
  
  
  # create results
  return(list(stat_res_pVal = result_sig_pVals, LFC_Values = LFC_table_per_lipidClass))
}
# Do PCA
doPCA <- function(dds,
                  colorTheme=NULL,
                  shapeVar,
                  colorVar,
                  xPC = "PC1",
                  yPC = "PC2"
                  ){
  #data <- as.data.frame(assay(vst(dds,blind=T)))
  data <- as.data.frame(assay(dds))
  pca <- prcomp(
    as.data.frame(t(data)),
    center = T,
    scale. = F
  )
  
  # how much variance is explained by each PC
  explVar <- pca$sdev^2/sum(pca$sdev^2)
  names(explVar) <- colnames(pca$x)
  # transform variance to percent
  percentVar <- round(100 * explVar, digits = 1)
  
  # Define data for plotting
  pcaData <- data.frame(pca$x,colData(dds))
  if(is.factor(pcaData[,colorVar])){
    #assume correct ordering of levels
  }else{
    pcaData[,colorVar] = factor(pcaData[,colorVar],
                                levels = as.character(unique(pcaData[,colorVar])),
                                ordered = T)
  }
if(is.null(colorTheme)){
  pca_plot <- ggplot(pcaData, aes(x = pcaData[,xPC],
                                  y = pcaData[,yPC],
                                  color=pcaData[,colorVar])) +
    geom_point(size =3,aes(shape=pcaData[,shapeVar]))+
    xlab(paste0(names(percentVar[xPC]),": ",percentVar[xPC], "% variance")) +
    ylab(paste0(names(percentVar[yPC]),": ", percentVar[yPC], "% variance")) +
    coord_fixed()+
    theme_classic()+
    theme(aspect.ratio = 1)
}else{
  pca_plot <- ggplot(pcaData, aes(x = pcaData[,xPC],
                                  y = pcaData[,yPC],
                                  color = pcaData[,colorVar]
                                  )) +
    geom_point(size =3,aes(shape=pcaData[,shapeVar]))+
    scale_color_manual(values = colorTheme,
                       name = colorVar)+
    xlab(paste0(names(percentVar[xPC]),": ",percentVar[xPC], "% variance")) +
    ylab(paste0(names(percentVar[yPC]),": ", percentVar[yPC], "% variance")) +
    coord_fixed()+
    theme_classic()+
    theme(aspect.ratio = 1)
}

  
  return(pca_plot)
}






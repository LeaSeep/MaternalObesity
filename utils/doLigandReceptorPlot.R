# Performs Ligand-Receptor Plot-Representation
# Input:
#  - take all DE Genes that are checked for being a ligand
#  - takes the respective ligand- dds Object (to get the annotation)
#  - takes all Genes that are checked for receptors
#  - data matrix from database that indicates ligand-receptor relationships

# Output:
#  - nothing
#  - saves plot directly to file


doLigandReceptorPlot <- function(
    allDEGenes_ligand,
    dds_obj_Ligand,
    allPresent_receptor,
    colorVar = "LFC",
    adjMatrix_LigandReceptor){
  
  cirocsData = adjMatrix_LigandReceptor[,c("ligand_ensembl_gene_id",
                                           "receptor_ensembl_gene_id",
                                           "ligand_gene_symbol",
                                           "receptor_gene_symbol")]
  # subset to provided data
  DE_ligand <- rownames(allDEGenes_ligand)[rownames(allDEGenes_ligand) %in% cirocsData$ligand_ensembl_gene_id]
  DE_receptor <- rownames(allPresent_receptor)[rownames(allPresent_receptor) %in% cirocsData$receptor_ensembl_gene_id]
  # check with connections are there
  cirocsData_present <- subset(cirocsData,
                               cirocsData$ligand_ensembl_gene_id %in% DE_ligand & 
                                 cirocsData$receptor_ensembl_gene_id %in%DE_receptor
                               )

  correspondingLFC_table = as.data.frame(allDEGenes_ligand[DE_ligand,colorVar])
  rownames(correspondingLFC_table) <- DE_ligand
  colnames(correspondingLFC_table) <- colorVar
  # add gene symbols
  correspondingLFC_table$SYMBOL <- rowData(dds_obj_Ligand)[rownames(correspondingLFC_table),"SYMBOL"]
  
  # Prep plot
  groupNames=c(
    rep("Kupffer Cell Ligands",length(cirocsData_present$ligand_gene_symbol)),
    rep("Hepatocyte receptor",length(cirocsData_present$receptor_gene_symbol))
    )
  names(groupNames)=c(cirocsData_present$ligand_gene_symbol,cirocsData_present$receptor_gene_symbol)
  
  grid.col=c(rep("#ed924c",length(cirocsData_present$ligand_gene_symbol)),
             rep("#754824",length(cirocsData_present$receptor_gene_symbol)))
  names(grid.col)=c(cirocsData_present$ligand_gene_symbol,cirocsData_present$receptor_gene_symbol)
  
  # Give more meaning to ordering & edge color
  cirocsData_present$LFC=0
  for(i in unique(cirocsData_present$ligand_gene_symbol)){
    cirocsData_present[which(cirocsData_present$ligand_gene_symbol==i),"LFC"] <- correspondingLFC_table[which(correspondingLFC_table$SYMBOL==i),colorVar]
  }
  
  
  cirocsData_present$LFC_noTrim <- cirocsData_present$LFC
  cirocsData_present[cirocsData_present$LFC>4,"LFC"]=4
  cirocsData_present[cirocsData_present$LFC<(-4),"LFC"]=-4
  
  col_fun = colorRamp2(
    c(max(cirocsData_present$LFC)*-1,
      -0.5,0,0.5,
      max(cirocsData_present$LFC)),
    c("darkred","bisque","white","azure","deepskyblue4")
    )
  
  cirocsData_present$linkColor=col_fun(cirocsData_present$LFC)
  
  OrderDF=c(
    cirocsData_present[order(cirocsData_present$LFC,decreasing = T),"ligand_gene_symbol"],
    cirocsData_present[order(cirocsData_present$LFC,decreasing = F),"receptor_gene_symbol"]
    )
  groupNames=groupNames[OrderDF]
  
  # do the plot
  pdf("LigandReceptor_Plot.pdf")
  circos.clear()
  
  chordDiagram(cirocsData_present[,c(3,4)],
               grid.col = grid.col, # kuppfer cells and hepatcytes colors
               transparency = 0.1,
               annotationTrack = c("grid"),
               grid.border=NA,
               preAllocateTracks = list(track.height = max(strwidth(unlist(c(cirocsData_present$ligand_gene_symbol,cirocsData_present$receptor_gene_symbol))))),
               group = groupNames,
               big.gap = 20,
               small.gap = 2,
               col=cirocsData_present$linkColor,
               direction.type = c("arrows"),
               directional = 1,
               link.arr.type = "big.arrow",
               link.border = NA,
               link.sort = "asis",
               target.prop.height = mm_h(3),
  )

  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  
  fine=seq(min(as.numeric(cirocsData_present$LFC)),-1*min(as.numeric(cirocsData_present$LFC)),0.2)
  fine_l=fine
  fine_l[1]=stringr::str_wrap("up in mat. obese (wt)",10)
  fine_l[length(fine_l)]=stringr::str_wrap("up in mat. lean (ko)",10)
  
  fine_l[c(2:10)]=""
  fine_l[c(12:20)]=""
  fine_l[c(22:30)]=""
  fine_l[c(32:40)]=""
  
  legend("topleft",
         inset=c(-0.025,0),
         bty = "n",
         legend=c(fine_l),
         fill=col_fun(fine),
         border=NA,
         y.intersp =-0.4,
         cex = 1)
  
  dev.off()
  
  return(cirocsData_present)
  
}

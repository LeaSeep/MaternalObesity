preprocessing <- function(dds_object, 
                          sumRowCountThres=10, 
                          protCodingOnly=T,
                          removeConstRows=T,
                          filterPerSample=F
                          ){
  data <- dds_object
  if(protCodingOnly){
    # take protein coding only if info is available
    if("GENETYPE" %in% colnames(rowData(data))){
      print("Subset to protein coding only")
      data <- data[rowData(data)$GENETYPE == "protein_coding",]
    }else{
      warning(
        "Did not find a colum 'GENETYPE' in 'rowData(dds_object)';
        no subsetting based on type done"
        )
    }
  }
  
  # remove anything constant
  if(removeConstRows){
    print("Remove all entities which are constant over all samples")
    data <- data[which(apply(as.data.frame(assay(data)),1,sd)!=0),]
  }
  
  if(filterPerSample){
    genes_to_keep <- rowSums(counts(data) >= sumRowCountThres) >= ceiling(0.25*ncol(data))
  }else{
    genes_to_keep <- rowSums(counts(data)) >= sumRowCountThres
  }
    

  dds_object_new <- data[genes_to_keep,]
  return(dds_object_new)
}
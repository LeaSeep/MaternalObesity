library(readxl)
library(SummarizedExperiment)

readInMetabolon <- function(
    filename,
    sheetname="Log Transformed Data",
    saveAsShinyRds=F,
    savingName=paste0("../data/",gsub("-","",Sys.Date()),"_Metabolomics")
    )
  {
  rawData=as.data.frame(readxl::read_excel(filename,
                                           trim_ws = T,
                                           sheet = sheetname))
  
  rownames(rawData)=rawData$PARENT_SAMPLE_NAME
  rawData$PARENT_SAMPLE_NAME=NULL
  colnames(rawData)=paste0("CHEMID_",colnames(rawData))
  
  #to get metablite x Samples
  rawData=as.data.frame(t(rawData))
  
  sampleData=as.data.frame(readxl::read_excel(filename,
                                              trim_ws = T,
                                              sheet = "Sample Meta Data"))
  
  rownames(sampleData)=sampleData$PARENT_SAMPLE_NAME
  sampleData$PARENT_SAMPLE_NAME=NULL
  
  rowAnno=as.data.frame(readxl::read_excel(filename,
                                           trim_ws = T,
                                           sheet = "Chemical Annotation"))
  
  rownames(rowAnno)=paste0("CHEMID_",rowAnno$CHEM_ID)
  
  object=list(type="Metabolomics",
              Matrix=rawData,
              sample_table=sampleData,
              annotation_rows=rowAnno)
  #rudimentary Checks
  all(rownames(rawData)==rownames(rowAnno))
  all(colnames(rawData)==rownames(sampleData))
  
  object_sumExp <- SummarizedExperiment(assays  = object$Matrix,
                                        rowData = object$annotation_rows[rownames(object$Matrix),],
                                        colData = object$sample_table)
  
  # check if correct and save as Shiny Object already?!
  if(saveAsShinyRds){
    toShiny=list()
    toShiny[["Metabolomics"]]=object
    toShiny[["Metabolomics_SumExp"]]=object_sumExp
    
    
    finalSaveNme=paste0(savingName,"_",sheetname,".rds")
    saveRDS(toShiny,finalSaveNme)
    print("Saved")
    return(toShiny)
  }else{
    return(object_sumExp)
  }
}

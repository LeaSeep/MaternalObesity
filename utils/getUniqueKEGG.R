getUnique <- function(DF, col ="KEGG"){
  allNAs <- which(!is.na(DF[,col]))
  print(
    paste0("Initially, ",nrow(DF[allNAs,]), " are not NA (from ",nrow(DF),")")
    )
  
  cleaned_first <- DF[allNAs,col]
  # remove all mets that cannot be uniquely adivsed
  multipleAssignment <- which(grepl(",",cleaned_first))
  print(
    paste0(
      "Another ",
      length(multipleAssignment),
      " are not uniquely, 1st one taken"
      )
    )
  
  for(i in multipleAssignment){
    cleaned_first[i] <- gsub(",.*","",cleaned_first[i])
  }
  
  return(cleaned_first)
}
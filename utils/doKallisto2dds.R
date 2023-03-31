doKallisto2dds <- function(
    kallisto_path,
    sampleTable,
    tx_anno_path
    ){

  tx_annotation <- read.delim(tx_anno_path, 
                              header = F , 
                              stringsAsFactors = F,
                              col.names = c("GENEID", "TXNAME", "SYMBOL", "GENETYPE"))
  
  # Define path where the Kallisto files are stored
  files <- paste(kallisto_path, sampleTable$ID, "/abundance.h5", sep = "")
  
  # Naming the entries in the vector assures correct column names in the expression tables
  names(files) <- rownames(sampleTable)
  
  # Import samples and perform the distribution of transcripts to genes
  tx_annotation <- tx_annotation[-1,]
  
  library(tximport)
  txi_kallisto <- tximport(files,
                           type="kallisto", 
                           tx2gene=tx_annotation[,2:1])
  
  library(DESeq2)
  
  dds_txi_hepatocytes <- DESeq2::DESeqDataSetFromTximport(txi = txi_kallisto, 
                                                          colData = sampleTable,
                                                          design = ~ Merged)
  # add gene annotation extracted from the gtf file
  gene_annotation <- tx_annotation[!duplicated(tx_annotation$GENEID),c("GENEID", "SYMBOL", "GENETYPE")]
  gene_annotation <- gene_annotation[match(rownames(dds_txi_hepatocytes), gene_annotation$GENEID),]
  
  rowData(dds_txi_hepatocytes) <- gene_annotation
  
  return(dds_txi_hepatocytes)
}
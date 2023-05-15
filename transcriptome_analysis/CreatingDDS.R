# Creating DDS object based on kallisto files
# can be downloaded from                         GEO ACCESSSION !!!!
# the results are saved in the data/ folder
source("../utils/doKallisto2dds.R")

if(gsub("^.+/","",getwd())!="transcriptome_analysis"){
  setwd("transcriptome_analysis")
}

# for both
dir <- "/Volumes/lseep@uni-bonn.de/SFB1454-P08/Mass_lab_maternal_obesity/Hif1a/RNA seq/RNA_seq/Data/output/kallisto/kallisto/"
tx_anno_path <- "../data/ID2SYMBOL_gencode_vM16_transcript.txt"

# HC ----
sample_table <- read.csv("../data/sampleTable_Hepatocytes_RNAseq.csv",row.names = 1)

dds_txi_hepatocytes <- doKallisto2dds(
  dir,
  sample_table,
  tx_anno_path
)

saveRDS(dds_txi_hepatocytes,"../data/DESeq_Obj_HC.rds")

# Save 3 Files for easy data-access if unfamiliar with DESeq Object
write.csv(
  as.data.frame(rowData(dds_txi_hepatocytes)),
  file = "../data/RowAnno_Hepatocytes_RNASeq_kallistoImport.csv"
  )
write.csv(
  as.data.frame(assay(dds_txi_hepatocytes)),
  file = "../data/CountMatrix_Hepatocytes_RNASeq_kallistoImport.csv"
)

# KC  ----
## Hif1a -----
sample_table <- read.csv("../data/sampleTable_KupfferCells_RNAseq.csv",row.names = 1)
dds_txi_KC <- doKallisto2dds(
  dir,
  sample_table,
  tx_anno_path
)

saveRDS(dds_txi_KC,"../data/DESeq_Obj_KC.rds")

# Save 3 Files for easy data-access if unfamiliar with DESeq Object
write.csv(
  as.data.frame(rowData(dds_txi_KC)),
  file = "../data/RowAnno_KupfferCells_RNASeq_kallistoImport.csv"
)
write.csv(
  as.data.frame(assay(dds_txi_KC)),
  file = "../data/CountMatrix_KupfferCells_RNASeq_kallistoImport.csv"
)

## WT ----
dir <- "/Volumes/lseep@uni-bonn.de/Maternal_obesity_paper_1_NB/WT/Adult_P76/RNA_seq/RNA_SMART_Mass_S160_C32/kallisto/kallisto/"
sample_table <- read.csv("../data/sampleTable_KupfferCells_WT_RNAseq.csv",row.names = 1)
sample_table <- subset(sample_table, subset = Age =="Adult")
sample_table$ID <- rownames(sample_table)

dds_txi_KC <- doKallisto2dds(
  dir,
  sample_table,
  tx_anno_path
)

saveRDS(dds_txi_KC,"../data/DESeq_Obj_KC_WT.rds")

# Save 3 Files for easy data-access if unfamiliar with DESeq Object
write.csv(
  as.data.frame(rowData(dds_txi_KC)),
  file = "../data/RowAnno_KupfferCells_WT_RNASeq_kallistoImport.csv"
)
write.csv(
  as.data.frame(assay(dds_txi_KC)),
  file = "../data/CountMatrix_KupfferCells_WT_RNASeq_kallistoImport.csv"
)

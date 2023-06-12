# Metabolome Analysis ----
## Overview ----
# This script includes the metabolome analysis of the serum of wt mice fed with
# different schemes (see manuscript). 
# Including Batch correction, statistical analysis and prep for Metaboanalyst.

## Input ----
# (all in data Folder):
#   - Metabolon output (Excel Workbook with multiple sheets)


## Output----
# (output to script's folder):
#   - ANOVA results for desired comparisons
#   - Visualisation of Batch corrected Data
#   - Metabolites to Metabridge

# SetUp ----
setwd("metabolome_analysis/")
source("../utils/readMetabolon.R")
source("../utils/doComBat.R")
source("../utils/doPCA.R")
source("../utils/getUniqueKEGG.R")
source("../utils/doANOVA_contrast.R")

output_result_list <- list()

colorTheme <- c("#c6c6c6","#606060","#c12c38","#e0775f","#f3b694","#fce2d0")

filename <- "../data/Metabolon_UHBO-06-21MD_DATA TABLES_extMetadata.xlsx"
checkObj <- readInMetabolon(
  filename,
  sheetname = "Log Transformed Data",
  saveAsShinyRds = F
  )
output_result_list[["Metabolomic_SumExp"]] <- checkObj
# ANOVA ----

## create Sets for Metaboanalyst ----
### HFD CD CD vs CD CD CD unadj. pVals ----

ANOVA_contrast_HFDCDCD_CDCDCD <- doANovaContrast(
  checkObj,
  design_fact1 = "MOTHER_DIET",
  design_fact2 = "OFFSPRING_DIET",
  contrast_nom = "HFD_CD:CD",
  contrast_ref = "CD_CD:CD"
  )
output_result_list[["ANOVA_HFDCDCD_vs_CDCDCD"]] <- ANOVA_contrast_HFDCDCD_CDCDCD

### Metaboanalyst Prep ----
# This are all metabolites where KEGG Id was provided 
# all Metabolites for which multiple IDs the first one is selected
# Metaboanalyst requests KEGGs as universe 

UniverseDF <- as.data.frame(rowData(checkObj))
UniverseKEGG <- getUnique(UniverseDF,col="KEGG")

write(UniverseKEGG,file = "KEGG_universe_for_Metabolyst.csv")

## create LookUP data set ----
# Get all significant Metabolites fo specified contrast (unadj. pVals)
# get HMDB (less missing values)
sigDF <- 
  ANOVA_contrast_HFDCDCD_CDCDCD[ANOVA_contrast_HFDCDCD_CDCDCD$`pVal_HFD_CD:CD_vs_CD_CD:CD`<0.05,]
sigDF <- as.data.frame(rowData(checkObj)[rownames(sigDF),])
sigKEGG <- getUnique(sigDF,col="KEGG")

write(sigKEGG,file = "KEGG_sigMetabolites_MaternalDiet_for_Metabolyst.csv")


# Batch Correction for DOE (day of experiment) ----
library(sva)
correctedObj <- doBatchCorrection(
  SumExp_obj = checkObj,
  design_factor = "GROUP_NAME",
  batch_factor = "day_of_experiment_day_when_serum_was_isolated"
  )

# PCA ----

pca_plot <- doPCA(correctedObj,
      colorTheme,
      "maternal_diet",
      "GROUP_NAME",
      xPC = "PC1",
      yPC = "PC2"
)

ggsave(filename = paste0("PCA_",Sys.Date(),".png"), plot=pca_plot)
ggsave(filename = paste0("PCA_",Sys.Date(),".svg"), plot=pca_plot)

# Save everything ----
saveRDS(output_result_list, file = "Metabolomics_results.rds")

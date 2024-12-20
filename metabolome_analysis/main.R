# Metabolome Analysis ----
## Overview ----
# This script includes the metabolome analysis of the serum of wt mice fed with
# different schemes (see manuscript). 
# Including Batch correction, statistical analysis and prep for Metaboanalyst.
# Overview of analysis can be found on the right (when in RStudio)

## Input ----
# (all in data Folder):
#   - Metabolon output (Excel Workbook with multiple sheets)


## Output----
# (output to script's folder):
#   - ANOVA results for desired comparisons
#   - Visualisation of Batch corrected Data
#   - PCA analysis


# SetUp ----
setwd("metabolome_analysis/")
source("../utils/readMetabolon.R")
source("../utils/doComBat.R")
source("../utils/doPCA.R")
source("../utils/getUniqueKEGG.R")
source("../utils/doANOVA_contrast.R")

output_result_list <- list()

colorTheme <- c("#c6c6c6","#606060","#c32b38","#fee2d1","#8a4094","#ebdeec")

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
correctedObj$GROUP_NAME <- factor(correctedObj$GROUP_NAME,
                                  levels = c("CD CD CD",
                                             "CD CD HFD",
                                             "HFD CD CD",
                                             "HFD HFD CD",
                                             "HFD CD HFD",
                                             "HFD HFD HFD"),
                                  ordered = T)

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

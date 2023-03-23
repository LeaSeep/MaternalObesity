# Transcriptomics Analysis ----
## Overview ----
# This script includes the transcritpome analysis of Hepatocytes (HC) as well as
# Kupffer-Cells (KC) - both Hif1a ko.
# Including DE, ORA, Ligand-Receptor and Co-expression analysis.

## Input ----
# (all in data Folder):
#   - Raw Count Matrix; produced as descriped in the manuscript, can also be
#     downloaded from GEO_ACCESSION
#   - sample annotation as provided, can also be downloaded from GEO_ACCESSION
#   - row annotation as provided, produced with `biomaRt::getBM` using ENSEMBL
#     Mus Musculus - downloaded on 05.07.2022 
#   - ligand-receptor interaction matrix (from CellTalk), parsed into rds object

## Output----
# (output to script's folder):
#   - DE-results object (KF & HC) - displayed in table-format in provided html file
#   - PCA plots
#   - Over-representation (ORA) plots for desired contrasts
#   - ORA-result objects
#   - Ligand-Receptor plot
#   - CoCena Heatmap
#   - CoCena Cluster ORA plots

# DE-Analysis ----

## KC ----

### PCA ----

## HC ----

### PCA ----

# ORA- Analysis ----

# Ligand-Receptor Analysis ----

# CoCena ----
## Co-ecpression analysis ----

## cluster ORA ----


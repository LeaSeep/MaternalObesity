# Project Description
This is the repository to redo all analysis done within the Maternal Obesity project.
(Manuscript Link will be added as soon as possible).
If you want to:
- redo the analysis for one omic type, please go to the respective subdirectory and execute the `main.R` script 
- check results tables regarding each Omic (DE analysis, ORA-analysis, CoCena) please checkout `Database.html`
- reuse custom functions, checkout the `utils`- directory

# Usage
## Run main analysis
To run the main analysis and produce all plots shown in the manuscript and supplementary you need to:
1. Open R
2. Set the working directory as the Folder that includes all analysis
3. Run `main.R` (manually or source it)

Note: If you source the file - current settings within the file might limit the output to the requested option.
To change you need to change respective settings (Marked in the scripts)

Example Lipidomics:
```{r}
setwd("MaternalObesity")
source("lipidome_analysis/main.R")
```
This will output the Heatmap of all Lipid-Classes for wt and the subset for ko.

## CoCena
The CoCena was conducted following closely the provided example Notebook. Therefore, to rerun the analysis one has to open the respective `.Rmd` Files.
In here, the code chunks can be either all executed at once or, again, manually.
To execute all, in RStudio, there should appear in the top left corner of the Editor a button "Run". Select here "Run All".

# Data Overview
## Data Repositories
The raw data with their respective Metadata can be found under the following addresses:
- Transcriptomics: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237408
- Lipidomics: http://dx.doi.org/10.21228/M81D9R
- Metabolomics: http://dx.doi.org/10.21228/M81D9R

## Processed Data
All data needed to run the analysis can be found in the `data` directory. Note, this might not necassarily be the raw data.
For transcriptomics the fastq files were summarized to a `.dds` object using the provided script `transcriptomics_analysis/CreatingDDS.R`.
For Lipidomics and Metabolomics, data was respectively parsed to three data matrices including the entite x samples matrices as well as one for sample Metadata and one for entitie Metadata. This format can be easily parsed to a `SummarizedExperiment`- object, which again can be reused in a plethora of Bioconductor Workflows.

## Database File
If you are just interested in more analysis-results, for example in exact statistical figures for a gene of interest you can rely on the `Database.html` file to check and also easily download respective data.
An html file can be opened in any browser. Please be patient, it will take a couple of minutes to load the file completely. On the left of the file you can find a dynamic look up table as soon as you start scrolling to jump to respective section. In this File you can find the data-tables concerning:
- Transcriptomics
  - DE analysis of Kupffer Cells and respective ORA analysis of the DE Sets (Up and Down regulated)
  - hCoCena identified cluster-gene sets and their respective ORA analysis ofr Kupffer Cells, Hepatocytes wt and ko
- Metabolomics ANOVA results for all metabolites
- Lipidomics 
  - Log Fold Changes actual numbers (no scaling, as done for visualisation purposes) for wt and ko
  - their respective actual p-values for wt and ko

# Further Information
If you need any help or have further questions about this repository please drop me an email!
lea.seep@uni-bonn.de or leaseep@gmail.com

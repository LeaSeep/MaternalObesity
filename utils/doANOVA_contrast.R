doANovaContrast <- function(
    SumExpObj,
    design_fact1 = "MOTHER_DIET",
    design_fact2 = "OFFSPRING_DIET",
    contrast_nom = "HFD_CD:CD", # Must be a present combi from levels of design_fact1 : design_fact2
    contrast_ref = "CD_CD:CD"
    ){
  library(contrast)
  library(data.table)
  dataMat <- as.data.frame(t(as.data.frame(assay(SumExpObj))))
  dataMat_merged=cbind(dataMat,colData(SumExpObj)[rownames(dataMat),c(design_fact1,design_fact2,"CLIENT_IDENTIFIER")])
  dataMat_merged[,design_fact1]=as.factor(dataMat_merged[,design_fact1])
  dataMat_merged[,design_fact2]=as.factor(dataMat_merged[,design_fact2])
  
  long <- melt(
    setDT(dataMat_merged),
    id.vars = c(design_fact1,design_fact2,"CLIENT_IDENTIFIER"), 
    variable.name = "metabolite",
    )
  long <- as.data.frame(long)
  long[,design_fact1]=as.factor(long[,design_fact1])
  long[,design_fact2]=as.factor(long[,design_fact2])
  
  # set Up recording df
  long[,paste0("pVal_",design_fact1)]=0
  long[,paste0("pVal_",design_fact2)]=0
  long$pVal_Int=0
  long[,paste0("pVal_",contrast_nom,"_vs_",contrast_ref)]=0
  
  design_formula <- as.formula(paste0("value~",design_fact1,"+",design_fact2,"+",design_fact1,":",design_fact2))
  contrast_fact1_part1=gsub(":.*","",contrast_nom)
  contrast_fact1_part2=gsub(".*:","",contrast_nom)
  contrast_fact2_part1=gsub(":.*","",contrast_ref)
  contrast_fact2_part2=gsub(".*:","",contrast_ref)
  
  for(i in unique(long$metabolite)){
    #anova done on loa data!!
    ANOVA_res=aov(formula = design_formula,data=subset(long,metabolite==i))
    long[long$metabolite==i,paste0("pVal_",design_fact1)]=unlist(summary(ANOVA_res))["Pr(>F)1"]
    long[long$metabolite==i,paste0("pVal_",design_fact2)]=unlist(summary(ANOVA_res))["Pr(>F)2"]
    long[long$metabolite==i,"pVal_Int"]=unlist(summary(ANOVA_res))["Pr(>F)3"]
    
    ANOVA_res=lm(formula = design_formula,data=subset(long,metabolite==i))
    listA <- list()
    listA[[design_fact1]] <- contrast_fact1_part1
    listA[[design_fact2]] <- contrast_fact1_part2
    
    listB <- list()
    listB[[design_fact1]] <- contrast_fact2_part1
    listB[[design_fact2]] <- contrast_fact2_part2
    
    contrast_test=contrast(ANOVA_res,
                           listA,
                           listB
                           )
    
    long[long$metabolite==i,paste0("pVal_",contrast_nom,"_vs_",contrast_ref)]=contrast_test$Pvalue
  }
  
  wide_df <- reshape(data=long,
                     idvar="metabolite",
                     v.names = "value",
                     timevar = c("CLIENT_IDENTIFIER"),
                     direction="wide")
  
  colnames(wide_df)=gsub("value.","",colnames(wide_df))
  rownames(wide_df) = wide_df$metabolite
  
  # add Biochmeical name
  wide_df = cbind(rowData(SumExpObj)[wide_df$metabolite,c("PATHWAY_SORTORDER","PLOT_NAME")],wide_df)
  
  wide_df[,paste0("pADJ_",contrast_nom,"_vs_",contrast_ref)]= p.adjust(wide_df[,paste0("pVal_",contrast_nom,"_vs_",contrast_ref)], method = "fdr") 
  
  return(wide_df)
}

# Batch Correction based on given design and specified variable as batch
# Input:
#   - Summarized Experiment object
#   - design formula
#   - name of factor which specifies the batch
# Output:
#   - batch-corrected Summarized Experiment object

doBatchCorrection <- function(
    SumExp_obj,
    design_factor,
    batch_factor){
  
  modcombat = model.matrix(as.formula(paste0("~",design_factor)) , data=colData(SumExp_obj))
  batch=colData(SumExp_obj)[,batch_factor]
  
  combat_edata = ComBat(
    dat=as.data.frame(assay(SumExp_obj)),
    batch=batch, 
    mod=modcombat, 
    par.prior=T, 
    prior.plots=F)
  
  assay(SumExp_obj) <- as.data.frame(combat_edata)
  
  return(SumExp_obj)
}
doOra <- function(
    geneSet,
    type=NULL,
    levelGOTerms=6,
    universe_entrez=NULL,
    filename){
  
  result = list()

  # assumes always ENSEMBL IDs
  geneSetChoice_tranlsated <- bitr(geneSet,
                                   fromType="ENSEMBL",
                                   toType="ENTREZID",
                                   OrgDb="org.Mm.eg.db")$ENTREZID
  if("GO" %in% type){
    EnrichmentRes_GO <- clusterProfiler::enrichGO(gene = geneSetChoice_tranlsated,
                                                  OrgDb="org.Mm.eg.db",
                                                  pvalueCutoff = 0.05,
                                                  ont = "BP",
                                                  universe = universe_entrez,
                                                  qvalueCutoff = 0.1,
                                                  readable=T)
    tryCatch(
      {
        EnrichmentRes_GO_filter = gofilter(EnrichmentRes_GO, level = levelGOTerms)
        EnrichmentRes_GO_filter_simple <- clusterProfiler::simplify(EnrichmentRes_GO_filter,
                                                                    cutoff = 0.8,
                                                                    by = "p.adjust",
                                                                    select_fun=min)
        },
      error=function(e){
        print("Simplyfing Go w.r.t to levels did not work, try to lower the number; the original object is returned and further processed")
        EnrichmentRes_GO_filter_simple <- EnrichmentRes_GO
      }
    )

    if(all(EnrichmentRes_GO_filter_simple@result$p.adjust>0.1)){
      print("GO nothing enriched - check result object")
    }else{
      png(paste0(filename,"_ORA","_GO.png"))
      print(dotplot(EnrichmentRes_GO_filter_simple))
      dev.off()
    }
    
    result[["GO"]] <- EnrichmentRes_GO_filter_simple@result

  }
  
  if("KEGG" %in% type){
    KEGGset <- msigdbr(
      species ="Mus musculus",
      category = "C2",
      subcategory = "KEGG"
    ) %>% dplyr::select(gs_name, entrez_gene)
    EnrichmentRes_Kegg <- clusterProfiler::enricher(
      gene = geneSetChoice_tranlsated,
      pvalueCutoff = 0.05,
      pAdjustMethod = "fdr",
      universe = universe_entrez,
      TERM2GENE = KEGGset
    )
    
    if(!is.null(EnrichmentRes_Kegg)){
      if(all(EnrichmentRes_Kegg@result$p.adjust>0.1)){
        print("Kegg nothing enriched - check result object")
      }else{
        png(paste0(filename,"_ORA","_KEGG.png"))
        print(dotplot(EnrichmentRes_Kegg))
        dev.off()
      }
      
      result[["KEGG"]] <- EnrichmentRes_Kegg@result
    }else{
      print("consider changing 'minGSSize' within function")
    }

  }
  
  if("HALLMARK" %in% type){
    Hallmarkset <- msigdbr(
      species ="Mus musculus",
      category = "H",
    ) %>% dplyr::select(gs_name, entrez_gene)
    EnrichmentRes_Hallmarks <- clusterProfiler::enricher(
      gene = geneSetChoice_tranlsated,
      pvalueCutoff = 0.05,
      pAdjustMethod = "fdr",
      universe = universe_entrez,
      TERM2GENE = Hallmarkset
    )
    
    if(!is.null(EnrichmentRes_Kegg)){    
      if(all(EnrichmentRes_Hallmarks@result$p.adjust>0.1)){
        print("Hallmark nothing enriched - check result object")
      }else{
        png(paste0(filename,"_ORA","_HALLMARK.png"))
        print(dotplot(EnrichmentRes_Hallmarks))
        dev.off()
    }
      result[["HALLMARK"]] <- EnrichmentRes_Hallmarks@result
    }else{
      print("consider changing 'minGSSize' within function")
    }
  }
  
  return(result)
}

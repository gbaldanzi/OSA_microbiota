# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlations with gut microbiota species using the extended model

  # This analysis only includes the species that were FDR significant in the main model 
  # Import species names identified in the main model including BMI
  #mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))
  

# Correlations

  res.extended.model <- lapply( exposures , spearman.function, 
                          x1 = outcomes,
                          covari = extended.model,
                          data = pheno)

  res.extended.model  <- do.call(rbind,res.extended.model)
  
  setDT(res.extended.model)
  setnames(res.extended.model,"x","MGS")
  

  # Save results for the extended model
  fwrite(res.extended.model, file = paste0(results.folder,"cor2_all.var_mgs.tsv"))
  
  
  # Species associated with OSA parameters after adjustment for BMI
  mgs.m1 <- unique(res.extended.model[q.value<.05,][["MGS"]]) 
  
  saveRDS(mgs.m1, paste0(results.folder,'mgs.m1.rds'))
  
  
  
  
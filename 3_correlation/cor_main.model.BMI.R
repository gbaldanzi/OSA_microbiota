# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlations with all main model covariates, including BMI
  
  res.main.model <- lapply(exposures,spearman.function, 
                x1=outcomes,
                covari = main.model.BMI,
                data = pheno)
  
  res.main.model <- do.call(rbind,res.main.model)
  
  setDT(res.main.model)
  setnames(res.main.model,"x","MGS")
  

  # Saving results 
  fwrite(res.main.model, file = paste0(results.folder,"cor.bmi_all.var_mgs.tsv"))
  
  # Species associated with OSA parameters after adjustment for BMI
  mgs.m1 <- unique(res.main.model[q.value<.05,][["MGS"]]) 
    
  saveRDS(mgs.m1, paste0(results.folder,'mgs.m1.rds'))
  
  
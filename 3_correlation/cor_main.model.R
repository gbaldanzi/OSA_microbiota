# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# Correlations using main model covariates except BMI
  res.main.model <- lapply(exposures,spearman.function, 
                x1=outcomes,
                covari = main.model,
                data = pheno)
  
  res.main.model <- do.call(rbind,res.main.model)
  
  setDT(res.main.model)
  setnames(res.main.model,"x","MGS")
  

  # Saving results 
  fwrite(res.main.model, file = paste0(results.folder,"cor_all.var_mgs.tsv"))
  
  mgs.fdr.main <- res.main.model[q.value<0.05, .(MGS,exposure)]

  saveRDS(mgs.fdr.main, paste0(results.folder, "mgs.fdr.mainmodel.rds"))
  
  
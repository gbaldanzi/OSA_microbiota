# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Import findings from main model w/o BMI

  mgs.fdr.main <- readRDS(paste0(results.folder, "mgs.fdr.mainmodel.rds"))

  mgs.fdr.main <- unique(mgs.fdr.main$MGS)
  
# Correlations with all main model covariates, including BMI

  
  res.main.model <- lapply(exposures,spearman.function, 
                x1 = mgs.fdr.main,
                covari = main.model.BMI,
                data = pheno)
  
  res.main.model <- do.call(rbind,res.main.model)
  
  setDT(res.main.model)
  setnames(res.main.model,"x","MGS")
  

  # Saving results 
  fwrite(res.main.model, file = paste0(results.folder,"cor.bmi_all.var_mgs.tsv"))
  
  
  
  
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlations with gut microbiota species using the "medication" model

  # This analysis only includes the species that were FDR significant in the main model 
  # Import species names identified in the main model including BMI
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))

# Correlations

  res.med.model <- lapply(exposures, spearman.function, 
                          x1=mgs.fdr,
                          covari = medication.model,
                          data = pheno)

  res.med.model  <- do.call(rbind,res.med.model)
  
  setDT(res.med.model)
  setnames(res.med.model,"x","MGS")
  

  
  # Save results for the full model
  fwrite(res.med.model, file = paste0(results.folder,"cor.med_all.var_mgs.tsv"))
  
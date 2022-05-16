# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# This script performs the sensitivity analysis by additional adjustment of medication use to the 
# extended model

# Correlations with gut microbiota species using the "medication" model

  # This analysis only includes the species that were FDR significant in the extended model 
  # Import species names identified in the extended model 
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))

# Correlations

  res.med.model <- lapply(exposures, spearman.function, 
                          x1 = mgs.fdr,
                          covari = medication.model,
                          data = pheno)

  res.med.model  <- do.call(rbind,res.med.model)
  
  setDT(res.med.model)
  setnames(res.med.model,"x","MGS")
  

  
  # Save results 
  fwrite(res.med.model, file = paste0(results.folder,"cor.med_all.var_mgs.tsv"))
  
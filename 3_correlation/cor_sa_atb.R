# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# This script performs the sensitivity analysis comparing the results 
# after exclusion of participants who have used antibiotic 

# This analysis only includes the species that were FDR significant in the extended model 

  # Remove participants that have used antibiotic in the last 6 months
  pheno.noatb <- pheno[atb6m=="no",]
  
  # Import species names identified in the extended model 
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))
  
  
  # Correlations 
  
  # This sensitivity analysis is adjusted for the extended model covariates

  res <- lapply(exposures,spearman.function, 
                          x1 = mgs.fdr,
                          covari = extended.model,
                          data = pheno.noatb)
  
  res  <- do.call(rbind,res)
  res$model= "sa_atb6m"
  
  setDT(res)
  setnames(res,"x","MGS")
  
  
  #Save results
  fwrite(res, file = paste0(results.folder,"corsaatb_all.var_mgs.tsv"))



# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script performs the sensitivity analysis comparing the results 
# after exclusion of participants who have self-reported CPOD, emphysema or chronic 
# bronchitis 

# This analysis only includes the species that were FDR significant in the extended model 

  # Import species names identified in the extended model 
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))

  # Remove participants that have used antibiotic in the last 6 months
  pheno.nolungdis <- pheno[-which(lungdisease=="yes"),]
  
  # Correlations 
  
  # This sensitivity analysis is adjusted for the extended model covariates

  res <- lapply(exposures,spearman.function, 
                          x1=mgs.fdr,
                          covari = extended.model,
                          data = pheno.nolungdis)
  
  res  <- do.call(rbind,res)
  res$model= "sa_lung"
  
  setDT(res)
  setnames(res,"x","MGS")
  
  
  #Save results
  fwrite(res, file = paste0(results.folder,"corsalung_all.var_mgs.tsv"))



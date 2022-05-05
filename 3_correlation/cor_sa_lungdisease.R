# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script carries out the sensitivity analysis comparing the results 
# after exclusion of participants who have self-reported CPOD, emphysema or chronic 
# bronchitis 

  # In this script, we will only include the species that were FDR significant in the main model 

  # This sensitivity analysis is adjusted for all main model covariates


  # Import species names identified in the main model including BMI
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))

  # Remove participants that have used antibiotic in the last 6 months
  pheno.nolungdis <- pheno[-which(lungdisease=="yes"),]
  
  # Correlations 

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



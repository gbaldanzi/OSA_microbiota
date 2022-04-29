# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# This script performs the sensitivity analysis comparing the results 
# after exclusion of participants who have used antibiotic 

# This sensitivity analysis is adjusted for all main model covariates

# This analysis only includes the species that were FDR significant in the main model 

  # Remove participants that have used antibiotic in the last 6 months
  pheno.noatb <- pheno[atb6m=="no",]
  
  # Import species names identified in the main model including BMI
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))
  
  
  # Correlations 

  res <- lapply(exposures,spearman.function, 
                          x1=mgs.fdr,
                          covari = main.model.BMI,
                          data = pheno.noatb)
  
  res  <- do.call(rbind,res)
  res$model= "sa_atb6m"
  
  setDT(res)
  setnames(res,"x","MGS")
  
  
  #Save results
  fwrite(res, file = paste0(results.folder,"corsaatb_all.var_mgs.tsv"))



# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Created: 2021-10-21

# This script was created to carry out the sensitivity analysis comparing the results 
# after exclusion of participants who have used antibiotic 

# In this script, we will only include the mgs that were FDR significant in the full model 

  # Import data 
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))

  # Remove participants that have used antibiotic in the last 6 months
  pheno <- pheno[atb6m=="no",]
  
  # Import signatures species 
  
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))  # Signature MGSs 
  
  
  # Correlations 

  res <- lapply(exposures,spearman.function, 
                          x1=mgs.m1,
                          covari = main.model.BMI,
                          data = pheno)
  
  res  <- do.call(rbind,res)
  res$model= "sa_atb6m"
  
  setDT(res)
  setnames(res,"x","MGS")
  
  
  #Save results
  fwrite(res, file = paste0(results.folder,"corsaatb_all.var_mgs.tsv"))



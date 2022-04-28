# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script carries out the sensitivity analysis comparing the results 
# after exclusion of participants who have self-reported CPOD, emphysema or chronic 
# bronchitis 

# In this script, we will only include the mgs that were FDR significant in the main model 

  # Import data 
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))

  # Remove participants that have used antibiotic in the last 6 months
  pheno <- pheno[-which(lungdisease=="yes"),]
  
  # Import signatures species 
  
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))  # Signature MGSs 
  
  
  # Correlations 

  res <- lapply(exposures,spearman.function, 
                          x1=mgs.fdr,
                          covari = main.model.BMI,
                          data = pheno)
  
  res  <- do.call(rbind,res)
  res$model= "sa_lung"
  
  setDT(res)
  setnames(res,"x","MGS")
  
  
  #Save results
  fwrite(res, file = paste0(results.folder,"corsalung_all.var_mgs.tsv"))



# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Created: 2021-10-21

# This script was created to carry out the sensitivity analysis comparing the results 
# after exclusion of participants who have used antibiotic 

# In this script, we will only include the mgs that were FDR significant in the full model 

  # Import data 
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

  # Remove participants that have used antibiotic in the last 6 months
  pheno <- pheno[atb6m=="no",]
  
  # Import signatures species 
  
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m2.rds'))  # Signature MGSs 
  
  
  # Correlations 

  
    #T90
  mgs.m1 = mgs.fdr$t90
  
  res.t90 <- spearman.function(x1=mgs.m1,x2="t90",
                               covari = full.model,
                               data = pheno[atb6m=="no",])
  
  #ODI
  mgs.m1 = mgs.fdr$odi
  
  res.odi <- spearman.function(x1=mgs.m1,x2="odi",
                               covari = full.model,
                               data = pheno[atb6m=="no",])
  
  #Poll results
  res <- rbind(res.t90, res.odi)
  res$model= "sa_atb6m"
  
  setDT(res)
  setnames(res,"x","MGS")
  
  
  #Save results
  fwrite(res, file = paste0(results.folder,"corsaatb_all.var_mgs.tsv"))



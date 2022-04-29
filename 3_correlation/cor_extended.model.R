# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlations with gut microbiota species using the extended model

  # This analysis only includes the species that were FDR significant in the main model 
  # Import species names identified in the main model including BMI
  mgs.fdr  = readRDS(paste0(results.folder,'mgs.m1.rds'))
  

# Correlations

  res.extended.model <- lapply(c("t90","odi"),spearman.function, 
                          x1=mgs.fdr,
                          covari = extended.model,
                          data = pheno)

  res.extended.model  <- do.call(rbind,res.extended.model)
  
  setDT(res.extended.model)
  setnames(res.extended.model,"x","MGS")
  

  # Save results for the extended model
  fwrite(res.extended.model, file = paste0(results.folder,"cor2_all.var_mgs.tsv"))
  
  
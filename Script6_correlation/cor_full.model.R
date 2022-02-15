# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-02

# Correlations with gut microbiota species using the full model

  # Import MGS identified in model 1
  # mgs.m1  = readRDS(paste0(results.folder,'mgs.m1.rds'))
  

# Correlations

  res.full.model <- lapply(exposures,spearman.function, 
                          x1=mgs.m1,
                          covari = full.model,
                          data = pheno)

  res.full.model  <- do.call(rbind,res.full.model)
  
  setDT(res.full.model)
  setnames(res.full.model,"x","MGS")
  
  
  # Merge with taxonomic annotation 
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  setnames(taxonomy,"maintax_mgs","MGS")
  
  res.full.model <- merge(res.full.model, taxonomy, by="MGS", all.x=T)
  
  
  # Save results for the full model
  fwrite(res.full.model, file = paste0(results.folder,"cor2_all.var_mgs.tsv"))
  
  
  # Signature species 
  res.full.model$q.value[res.full.model$q.value>=0.001] <- round(res.full.model$q.value[res.full.model$q.value>=0.001] , digits = 3)
  
  setDT(res.full.model)
  
  mgs.fdr.ahi = res.full.model[exposure=="ahi" & q.value<.05 ,MGS]
  mgs.fdr.t90 = res.full.model[exposure=='t90' & q.value<.05 ,MGS]
  mgs.fdr.odi = res.full.model[exposure=='odi' & q.value<.05 ,MGS]
                     
  
  mgs.m2 <- list(mgs.fdr.ahi = mgs.fdr.ahi,
                 mgs.fdr.t90 = mgs.fdr.t90,
                 mgs.fdr.odi = mgs.fdr.odi) 
  
  
  saveRDS(mgs.m2, paste0(results.folder,'mgs.m2.rds'))
 
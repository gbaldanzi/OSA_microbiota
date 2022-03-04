# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-22

# Correlations with gut microbiota species using the extended model

  # Import MGS identified in model 1
  # mgs.m1  = readRDS(paste0(results.folder,'mgs.m1.rds'))
  

# Correlations

  res.extended.model <- lapply(exposures,spearman.function, 
                          x1=outcomes,
                          covari = extended.model,
                          data = pheno)

  res.extended.model  <- do.call(rbind,res.extended.model)
  
  setDT(res.extended.model)
  setnames(res.extended.model,"x","MGS")
  
  
  # Merge with taxonomic annotation 
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  setnames(taxonomy,"maintax_mgs","MGS")
  
  res.extended.model <- merge(res.extended.model, taxonomy, by="MGS", all.x=T)
  
  
  # Save results for the extended model
  fwrite(res.extended.model, file = paste0(results.folder,"cor2_all.var_mgs.tsv"))
  
  
  # Signature species 
  #setDT(res.extended.model)
  
  #mgs.fdr.ahi = res.extended.model[exposure=="ahi" & q.value<.05 ,MGS]
  #mgs.fdr.t90 = res.extended.model[exposure=='t90' & q.value<.05 ,MGS]
  #mgs.fdr.odi = res.extended.model[exposure=='odi' & q.value<.05 ,MGS]
                     
  
  #mgs.m2 <- list(mgs.fdr.ahi = mgs.fdr.ahi,
  #               mgs.fdr.t90 = mgs.fdr.t90,
  #               mgs.fdr.odi = mgs.fdr.odi) 
  
  
  #saveRDS(mgs.m2, paste0(results.folder,'mgs.m2.rds'))
 
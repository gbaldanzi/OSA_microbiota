# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-22

# Correlations 
  
  res.main.model <- lapply(exposures,spearman.function, 
                x1=outcomes,
                covari = main.model.BMI,
                data = pheno)
  
  res.main.model <- do.call(rbind,res.main.model)
  
  setDT(res.main.model)
  setnames(res.main.model,"x","MGS")
  
  # Merge with taxonomic annotation 
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  setnames(taxonomy,"maintax_mgs","MGS")
  
  res.main.model <- merge(res.main.model, taxonomy, by="MGS", all.x=T)

  # Saving results 
  fwrite(res.main.model, file = paste0(results.folder,"cor.bmi_all.var_mgs.tsv"))
  
  mgs.m1 <- unique(res.main.model[q.value<.05,][["MGS"]]) 
    
  saveRDS(mgs.m1, paste0(results.folder,'mgs.m1.rds'))
  
  
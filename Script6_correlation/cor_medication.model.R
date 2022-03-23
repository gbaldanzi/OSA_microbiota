# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-22

# Correlations with gut microbiota species using the "medication" model

  # Import MGS identified in the main mode
  mgs.m1  = readRDS(paste0(results.folder,'mgs.m1.rds'))
  

# Correlations

  res.med.model <- lapply(exposures,spearman.function, 
                          x1=mgs.m1,
                          covari = medication.model,
                          data = pheno)

  res.med.model  <- do.call(rbind,res.med.model)
  
  setDT(res.med.model)
  setnames(res.med.model,"x","MGS")
  
  
  # Merge with taxonomic annotation 
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  setnames(taxonomy,"maintax_mgs","MGS")
  
  res.med.model <- merge(res.med.model, taxonomy, by="MGS", all.x=T)
  
  
  # Save results for the full model
  fwrite(res.med.model, file = paste0(results.folder,"cor.med_all.var_mgs.tsv"))
  
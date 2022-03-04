# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-22

# Correlations 
  res.main.model <- lapply(exposures,spearman.function, 
                x1=outcomes,
                covari = main.model,
                data = pheno)
  
  res.main.model <- do.call(rbind,res.main.model)
  
  setDT(res.main.model)
  setnames(res.main.model,"x","MGS")
  
  # Merge with taxonomic annotation 
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  setnames(taxonomy,"maintax_mgs","MGS")
  
  res.main.model <- merge(res.main.model, taxonomy, by="MGS", all.x=T)

  # Saving results 
  fwrite(res.main.model, file = paste0(results.folder,"cor_all.var_mgs.tsv"))

  
  
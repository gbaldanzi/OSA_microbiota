# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


  # Import findings from main model w/o BMI

  mgs.fdr.main <- readRDS(paste0(results.folder, "mgs.fdr.mainmodel.rds"))

  mgs.fdr.main <- unique(mgs.fdr.main$MGS)

  # Correlations with gut microbiota species using the extended model

# Correlations

  res.extended.model <- lapply( exposures , spearman.function, 
                          x1 = outcomes,
                          covari = extended.model,
                          data = pheno)

  res.extended.model  <- do.call(rbind,res.extended.model)
  
  setDT(res.extended.model)
  setnames(res.extended.model,"x","MGS")
  
  
  # FDR-adjustment 
  res.extended.model[, q.value := NA]
  for(exp in c("ahi","t90","odi")) {
  res.extended.model[MGS %in% mgs.fdr.main & exposure == exp, q.value := p.adjust(p.value,method="BH")]
  }
  res.extended.model[q.value>=0.001, q.value := round(q.value,3)]
  

  # Save results for the extended model
  fwrite(res.extended.model, file = paste0(results.folder,"cor2_all.var_mgs.tsv"))
  
  
  # Species associated with OSA parameters after adjustment for BMI
  mgs.m1 <- unique(res.extended.model[q.value<.05,][["MGS"]]) 
  
  saveRDS(mgs.m1, paste0(results.folder,'mgs.m1.rds'))
  
  
  
  
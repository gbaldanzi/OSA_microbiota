# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-01

# Removing rare species (species present in less than 5 individuals)   
  # Calculating MGS prevalence 
  species.names=grep("____",names(pheno),value=T)
# presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
  data_pa <- decostand(x = pheno[,species.names,with=F], "pa")
# calculate sum per species
  data_sum <- data.frame(prevalence=apply(data_pa, 2, sum)/nrow(data_pa))
  data_sum$MGS = rownames(data_sum)
  
  a = data_sum$MGS[data_sum$prevalence <= 1/100] # 1,602 species have a prevalence greater than 1%
  pheno <- pheno[ , -a, with=F] 
  
  
# Correlations 
  
  outcomes  <-  grep("___",names(pheno),value=T)
  exposures <- c("ahi","t90","odi")
  
  res.basic.model <- lapply(exposures,spearman.function, 
                x1=outcomes,
                covari = basic.model,
                data = pheno)
  
  res.basic.model <- do.call(rbind,res.basic.model)
  
  setDT(res.basic.model)
  setnames(res.basic.model,"x","MGS")

  # Saving results 
  fwrite(res.basic.model, file = paste0(results.folder,"cor_all.var_mgs_temp.tsv"))
  
  mgs.m1 <- unique(res.basic.model[q.value<.05,][["MGS"]]) 
    
  saveRDS(mgs.m1, paste0(results.folder,'mgs.m1_temp.rds'))
  
  
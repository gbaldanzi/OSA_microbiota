# Project: Sleep apnea and gut microbiota 
# Gabriel Baldanzi

# Enrichment analysis for gut metabolic modules (GMM) using extended model ranked p-values stratified by 
# the direction of the Spearman's correlation coefficients

rm(list=ls())

  # Loading packages 
  library(data.table)
  library(tidyr)
  library(fgsea)
  library(stringr)
  
  
  # Functions 
  source("0_functions/MGS.Enrich.function.R")

  # Folders 
  input <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/data_processed/'
  results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
  work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'

  
  # Import Spearman's correlation results from extended model 
  res <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
  
  res[,mgs:=cutlast(MGS,9)]
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   ODI = res[exposure=="odi",])
  

  # List of modules ####
  load(paste0(input,'MGS_HG3A.GMMs2MGS.RData')) # object = MGS_HG3A.GMMs2MGS
  
  # Enrichment analysis   
  
  ## Positive correlations ####
  
  list.modules <-  MGS_HG3A.GMMs2MGS
  
  for(m in names(list.modules)){
    list.modules[[m]] <- list.modules[[m]][list.modules[[m]] %in% res[rho>0,mgs]]
  }

    res.pos = list()[1:3]
    
    for(i in 1:3){
    
      res.pos[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                       p.value.name="p.value",
                       MGS.var.name = "mgs",
                       direction = "positive",
                       enrich.var.list = list.modules)
      
      }
      
    res.pos <- do.call(rbind,res.pos)
    

      ## Negative correlations ####

      list.modules <-  MGS_HG3A.GMMs2MGS
        
      for(m in names(list.modules)){
          list.modules[[m]] <- list.modules[[m]][list.modules[[m]] %in% res[rho<0,mgs]]
      }
        
      res.neg = list()[1:3]
        
      for(i in 1:3){
          
        res.neg[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                           p.value.name="p.value",
                           MGS.var.name = "mgs",
                           enrich.var.list = list.modules,
                           direction = "negative")
      }
        
      res.neg <- do.call(rbind, res.neg)
        
      # Save results 
      res <- rbind(res.pos,res.neg)
        
      fwrite(res, file = paste0(results.folder,"ea_GMM.tsv"))
        
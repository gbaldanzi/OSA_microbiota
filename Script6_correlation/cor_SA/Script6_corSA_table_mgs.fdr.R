# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-09-14

# This code will create a table with the results for the SA for the relevant MGS

  library(data.table)
  #library(tidyverse)

  # Importing results 
  res.m2 <-  fread("/home/baldanzi/Sleep_apnea/Results/cor2_all.var_mgs.tsv")
  res.sa <-  fread("/home/baldanzi/Sleep_apnea/Results/corsa_all.var_mgs.tsv")
  
  res.sa[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  mgs.bmi <- res.m2[q.value<0.05 & exposure=="BMI",MGS]
  mgs.fdr <- unique(res.m2[q.value<0.05 & exposure %in% c("ahi","t90"),MGS])
  mgs.fdr <- mgs.fdr[!mgs.fdr %in% mgs.bmi]
  
  # Keep only the significant MGS
  a <- c("MGS","cor.coefficient", "p.value","q.value")
  res.saahi <- res.sa[MGS %in% mgs.fdr & exposure=="ahi", a, with = F]
  res.sat90 <- res.sa[MGS %in% mgs.fdr & exposure=="t90", a, with = F]
  a <- c("MGS","cor.coefficient", "p.value")
  res.m2ahi <- res.m2[MGS %in% mgs.fdr & exposure == "ahi", a, with = F]
  res.m2t90 <- res.m2[MGS %in% mgs.fdr & exposure == "t90", a, with = F]
  
  res.list <- list(m2ahi= res.m2ahi, saahi=res.saahi,
                   m2t90=res.m2t90, sat90=res.sat90)
  
  name.cols <- NULL
  
  res.list[[2]] <- res.list[[2]][q.value>.001, q.value:=round(q.value,3)]
  res.list[[4]] <- res.list[[4]][q.value>.001, q.value:=round(q.value,3)]
  
  for(i in 1:4){
    setnames(res.list[[i]], "cor.coefficient", "cor")
    
    res.list[[i]] <- res.list[[i]][order(MGS)]
    
    res.list[[i]][,MGS:=NULL]
    
    res.list[[i]] <- res.list[[i]][p.value>.001, p.value:=round(p.value,3)]
    res.list[[i]] <- res.list[[i]][, cor:=round(cor,3)]
    
    names(res.list[[i]]) <- paste0(names(res.list[[i]]),"_",names(res.list[i]))
    
    name.cols <- c(name.cols, names(res.list[[i]]))
    
  }
  

  
  res.table <- do.call(cbind,res.list)
  names(res.table) <- name.cols
  res.table$MGS <- mgs.fdr[order(mgs.fdr)]
  
  res.table <- res.table[,c(11,1:10)]
  
  saveRDS(res.table,file="/home/baldanzi/Sleep_apnea/Results/table.res_sa_mgs.fdr.rds")
  
  
  
  
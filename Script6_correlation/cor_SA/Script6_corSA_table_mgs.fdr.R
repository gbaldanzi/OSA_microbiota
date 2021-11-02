# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-09-14

# This code will create a table with the results for the SA for the relevant MGS

  library(data.table)
  library(tidyverse)

  # Importing results 
  res.m2 <-  fread("/home/baldanzi/Sleep_apnea/Results/cor2_all.var_mgs.tsv")
  res.sa <-  fread("/home/baldanzi/Sleep_apnea/Results/corsa_all.var_mgs.tsv")
  
  res.sa[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  res.sa[p.value>=0.001, p.value:=round(p.value, digits = 3)]
  res.m2[p.value>=0.001, p.value:=round(p.value, digits = 3)]
  
  res.sa[, rho:=round(rho, digits = 3)]
  res.m2[, rho:=round(rho, digits = 3)]
  
  res.sa[,rho:=cor.coefficient]
  res.m2[,rho:=cor.coefficient]
  
  mgs.fdr <- readRDS("/home/baldanzi/Sleep_apnea/Results/mgs.m2.rds")
  
  # Keep only the significant MGS
  a <- c("MGS","rho", "p.value","q.value","model")
  res.saahi <- res.sa[MGS %in% mgs.fdr$mgs.fdr.ahi & exposure=="ahi", a, with = F]
  res.sat90 <- res.sa[MGS %in% mgs.fdr$mgs.fdr.t90 & exposure=="t90", a, with = F]
  res.m2ahi <- res.m2[MGS %in% mgs.fdr$mgs.fdr.ahi & exposure == "ahi", a, with = F]
  res.m2t90 <- res.m2[MGS %in% mgs.fdr$mgs.fdr.t90 & exposure == "t90", a, with = F]
  
  res.ahi <- rbind(res.m2ahi, res.saahi) %>% 
    pivot_wider(id_cols = MGS, names_from=model, values_from = c("rho", "p.value" ,"q.value"))                                                     
  
  
  res.t90 <- rbind(res.m2t90, res.sat90) %>% 
    pivot_wider(id_cols = MGS, names_from=model, values_from = c(rho, p.value,q.value))                                                     
  
  list.res <- list(ahi = res.ahi, t90 = res.t90)

  
  saveRDS(list.res,file="/home/baldanzi/Sleep_apnea/Results/table.res_sa_mgs.fdr.rds")
  
  
  
  
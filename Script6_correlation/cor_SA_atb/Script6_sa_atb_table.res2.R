# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-10-21

# Lastest update: 2021-10-21

# This code will produce a table from the correlations of AHI, BMI, 
# and T90 with MGSs in the Sensitivity Analysis with Antibiotic 

pacman::p_load(data.table,tidyr,dplyr)

  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"

   # Results from sensitivity analysis 
  # Importing results 
  res.m2 <-  fread("/home/baldanzi/Sleep_apnea/Results/cor2_all.var_mgs.tsv")
  res.sa <-  fread(paste0(input,"cor_sa_atb6m_step2_all.var_mgs.tsv"))
  
  res.sa[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  res.sa[p.value>=0.001, p.value:=round(p.value, digits = 3)]
  res.m2[p.value>=0.001, p.value:=round(p.value, digits = 3)]
  
  res.m2[,rho:=cor.coefficient]
  
  res.sa[, rho:=round(rho, digits = 3)]
  res.m2[, rho:=round(rho, digits = 3)]
  
  
  
  mgs.fdr <- readRDS("/home/baldanzi/Sleep_apnea/Results/mgs.m2.rds")
  
  # Keep only the significant MGS
  a <- c("MGS","rho", "p.value","q.value","model","N")
  res.saahi <- res.sa[MGS %in% mgs.fdr$mgs.fdr.ahi & exposure=="ahi", a, with = F]
  res.sat90 <- res.sa[MGS %in% mgs.fdr$mgs.fdr.t90 & exposure=="t90", a, with = F]
  res.m2ahi <- res.m2[MGS %in% mgs.fdr$mgs.fdr.ahi & exposure == "ahi", a, with = F]
  res.m2t90 <- res.m2[MGS %in% mgs.fdr$mgs.fdr.t90 & exposure == "t90", a, with = F]
  
  res.ahi <- rbind(res.m2ahi, res.saahi) %>% 
    pivot_wider(id_cols = MGS, names_from=model, values_from = c(rho, p.value ,q.value, N))                                                     
  
  
  res.t90 <- rbind(res.m2t90, res.sat90) %>% 
    pivot_wider(id_cols = MGS, names_from=model, values_from = c(rho, p.value,q.value, N))  
  
  
  setDT(res.ahi)
  setDT(res.t90)
  
  list.res <- list(ahi = res.ahi, t90 = res.t90)
  
  
  # Save list of table  table 
  saveRDS(list.res,file="/home/baldanzi/Sleep_apnea/Results/table.res_sa_atb_step2.rds")
  
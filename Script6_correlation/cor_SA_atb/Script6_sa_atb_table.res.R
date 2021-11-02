# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-10-21

# Lastest update: 2021-10-21

# This code will produce a table from the correlations of AHI, BMI, 
# and T90 with MGSs in the Sensitivity Analysis with Antibiotic 

pacman::p_load(data.table,tidyr,dplyr)
source("Script6.Functions.R")

  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"

   # Results from sensitivity analysis 
  res <- fread(paste0(input,"cor_sa_atb3m_all.var_mgs.tsv"))
  
  res[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  mgs.fdr <- unique(res[q.value<0.05, MGS])
  
  res <- res[MGS %in% mgs.fdr,]
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   BMI = res[exposure=="BMI",])
  
  res.list = lapply(res.list,Clean.Correlation.Results)
  

  # From long to wide
  res.df = rbind(res.list[[1]],res.list[[2]],res.list[[3]])
  #a = c("cor.","p.value","q.value","N")
  #res.df = dcast(setDT(res.df), MGS~exposure,value.var=a)
  
  res.df <- res.df %>% pivot_wider(id_cols = MGS, names_from=exposure, 
                                   values_from= c("rho","p.value","q.value","N"))
  
  # Save table 
  saveRDS(res.df,file="/home/baldanzi/Sleep_apnea/Results/table.res_sa_atb.rds")
  
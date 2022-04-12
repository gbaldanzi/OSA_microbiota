# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-18

# Inferential Statistics 

# This code will produce a table from the correlations of AHI, BMI, 
# and T90 with MGSs in model 2

pacman::p_load(data.table,ggplot2,ggvenn, tidyr,dplyr)
source("Script6.Functions.R")

# input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"

# Import data

  # Results from model1 
  res <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   ODI = res[exposure=="odi",],
                   BMI = res[exposure=="BMI",])
  
  res.list = lapply(res.list,Clean.Correlation.Results)
  
  # From long to wide
  res.df2 = do.call(rbind,res.list)
 
  res.df2 <- res.df2 %>% pivot_wider(id_cols = MGS, names_from=exposure, 
                                   values_from= c("rho","p.value","q.value","N"))
  
  # Save table 
  saveRDS(res.df2,file="/home/baldanzi/Sleep_apnea/Results/table.res2.rds")
  

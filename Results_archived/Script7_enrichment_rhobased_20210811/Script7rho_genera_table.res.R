# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-05

# Inferential Statistics 

# This code will produce a table with overrepresentation analysis of genera
# correlated to AHI, T90 and BMI. -

# Based on rho


  pacman::p_load(data.table)
  
  # input and output folders 
  input <-  "/home/baldanzi/Sleep_apnea/Results/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"
  
  #Importing
  res <- readRDS(file = paste0(input,"ea_genera_rho.rds"))

  #Cleaning
  res <-  lapply(res,function(x){x[order(q.value) & q.value<0.05,.(genera,pval,q.value)]})

  #Saving 
  saveRDS(res,file=paste0(output,"ea_genera_table_rho.rds"))

  
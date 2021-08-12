# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-05

# Inferential Statistics 

# This code will produce a table with overrepresentation analysis of modules
# correlated to AHI, T90 and BMI


  pacman::p_load(data.table)
  
  # input and output folders 
  input <-  "/home/baldanzi/Sleep_apnea/Results/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"
  
  #Importing
  res.rho <- fread(paste0(input,"ea_modules_rho.tsv"))

  #Cleaning
  res.rho <-  res.rho[order(q.value) & q.value<0.05,.(exposure,type,category,name,pval,q.value)]
  
  #Saving 
  saveRDS(res.rho,file=paste0(output,"ea_modules_table_rho.rds"))

  
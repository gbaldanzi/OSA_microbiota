# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-05

# Inferential Statistics 

# This code will produce a table with overrepresentation analysis of subpathways
# among the correlated metabolites  for to AHI, T90 and BMI

  rm(list=ls())
  pacman::p_load(data.table)
  
  # input and output folders 
  input <-  "/home/baldanzi/Sleep_apnea/Results/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"
  
  #Importing
  res <- fread(paste0(input,"ea_subpath_rho.tsv"))

  #Cleaning
  res[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  res[, p.value:=ifelse(pval>=0.001,round(pval, digits = 3),pval)]
  res <-  res[order(q.value) & q.value<0.05,.(exposure,pathway,p.value,q.value)]

  
  #Saving 
  saveRDS(res,file=paste0(output,"ea_subpath_table_rho.rds"))

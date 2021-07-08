# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-05

# Inferential Statistics 

# This code will produce a table with overrepresentation analysis of family
# correlated to AHI, T90 and BMI

# Based on rho

  rm(list=ls())
  pacman::p_load(data.table)
  
  # input and output folders 
  input <-  "/home/baldanzi/Sleep_apnea/Results/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"
  
  #Importing
  res.both <- readRDS(file = paste0(input,"ea_family_rho.rds"))
  
  #Cleaning
  res.both <-  lapply(res.both,function(x){x[order(q.value) & q.value<0.05,.(family,pval,q.value)]})

  #Saving 
  saveRDS(res.both,file=paste0(output,"ea_family_table_rho.rds"))

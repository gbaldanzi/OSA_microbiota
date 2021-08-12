# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-01

# Inferential Statistics 

# This code will produce a table with overrepresentation analysis of genera
# correlated to AHI, T90 and BMI


  pacman::p_load(data.table)
  
  # input and output folders 
  input <-  "/home/baldanzi/Sleep_apnea/Results/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"
  
  #Importing
  res.both <- readRDS(file = paste0(input,"ea_genera_both.rds"))
  res.pos <- readRDS(file = paste0(input,"ea_genera_pos.rds"))
  res.neg <- readRDS(file = paste0(input,"ea_genera_neg.rds"))
  
  #Cleaning
  res.both <-  lapply(res.both,function(x){x[order(q.value) & q.value<0.05,.(genera,pval,q.value)]})
  res.pos <-  lapply(res.pos,function(x){x[order(q.value) & q.value<0.05,.(genera,pval,q.value)]})
  res.neg <-  lapply(res.neg,function(x){x[order(q.value) & q.value<0.05,.(genera,pval,q.value)]})
  
  #Saving 
  saveRDS(res.both,file=paste0(output,"ea_genera_table_both.rds"))
  saveRDS(res.pos,file=paste0(output,"ea_genera_table_pos.rds"))
  saveRDS(res.neg,file=paste0(output,"ea_genera_table_neg.rds"))

  
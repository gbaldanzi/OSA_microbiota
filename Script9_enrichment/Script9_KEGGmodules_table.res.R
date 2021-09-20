# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-01

# Inferential Statistics 

# This code will produce a table with overrepresentation analysis of modules
# correlated to AHI, T90 and BMI


  pacman::p_load(data.table)
  
  # input and output folders 
  input <-  "/home/baldanzi/Sleep_apnea/Results/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"
  
  #Importing
  res.both <- fread(paste0(input,"ea_modules_both.tsv"))
  res.pos <- fread(paste0(input,"ea_modules_pos.tsv"))
  res.neg <- fread(paste0(input,"ea_modules_neg.tsv"))
  
  #Cleaning
  res.both[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  res.pos[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  res.neg[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  res.both <-  res.both[order(q.value) & q.value<0.05,.(exposure,type,category,name,pval,q.value)]
  res.pos <-  res.pos[order(q.value) & q.value<0.05,.(exposure,type,category,name,pval,q.value)]
  res.neg <-  res.neg[order(q.value) & q.value<0.05,.(exposure,type,category,name,pval,q.value)]

  #Saving 
  saveRDS(res.both,file=paste0(output,"ea_modules_table_both.rds"))
  saveRDS(res.pos,file=paste0(output,"ea_modules_table_pos.rds"))
  saveRDS(res.neg,file=paste0(output,"ea_modules_table_neg.rds"))

  
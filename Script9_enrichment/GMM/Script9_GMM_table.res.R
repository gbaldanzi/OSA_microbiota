# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-10-30

# Inferential Statistics 

# This code will produce a table with overrepresentation analysis of GMM modules
# correlated to AHI, T90 and BMI


  pacman::p_load(data.table)
  
  # input and output folders 
  input <-  "/home/baldanzi/Sleep_apnea/Results/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"
  
  #Importing
  res.pos <- fread(paste0(input,"ea_GMM_pos.tsv"))
  res.neg <- fread(paste0(input,"ea_GMM_neg.tsv"))
  
  # Import GMM names 
  gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  
  res.pos <- merge(res.pos, gmm.names, by.x="pathway", by.y="Module")
  res.neg <- merge(res.neg, gmm.names, by.x="pathway", by.y="Module")
  
  #Cleaning
  res.pos[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  res.neg[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  res.pos[,NES:=round(NES,3)]
  res.neg[,NES:=round(NES,3)]
  
  res.pos <-  res.pos[order(q.value) & q.value<0.05,.(exposure,pathway,Name,HL1, HL2,NES,pval,q.value)]
  res.neg <-  res.neg[order(q.value) & q.value<0.05,.(exposure,pathway,Name,HL1, HL2,NES,pval,q.value)]

  #Saving 
  saveRDS(res.pos,file=paste0(output,"ea_GMM_table_pos.rds"))
  saveRDS(res.neg,file=paste0(output,"ea_GMM_table_neg.rds"))

  
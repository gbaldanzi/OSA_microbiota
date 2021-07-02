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

  # MGSs identified in model 2
  mgs.m2 <- readRDS('/home/baldanzi/Sleep_apnea/Results/mgs.m2.rds')
  
  # Results from model1 
  res.ahi <- fread(paste0(input,"cor2_ahi_mgs.tsv"))
  res.bmi <- fread(paste0(input,"cor2_BMI_mgs.tsv"))
  res.t90 <- fread(paste0(input,"cor2_t90_mgs.tsv"))
  
  res.list = list(res.ahi,res.bmi,res.t90)
  res.list = lapply(res.list,Clean.Correlation.Results)
  
  # From long to wide
  res.df2 = rbind(res.list[[1]],res.list[[2]],res.list[[3]])
  a = c("cor.","p.value","q.value","N")
  res.df2 = dcast(setDT(res.df2), MGS~exposure,value.var=a)
  
  # Save table 
  saveRDS(res.df2,file="/home/baldanzi/Sleep_apnea/Results/table.res2.rds")
  

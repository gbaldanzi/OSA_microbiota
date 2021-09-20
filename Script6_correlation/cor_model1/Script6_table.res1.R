# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-17

# Lastest update: 2021-08-02

# Inferential Statistics 

# This code will produce a table from the correlations of AHI, BMI, 
# and T90 with MGSs in model 1

pacman::p_load(data.table,ggplot2,tidyr,dplyr)
source("Script6.Functions.R")

# input and output folders 
input = "/home/baldanzi/Sleep_apnea/Results/"

# Import data

  # MGSs identified in model 1
  mgs.m1 <- readRDS('/home/baldanzi/Sleep_apnea/Results/mgs.m1.rds')
  
  # Results from model1 
  res <- fread(paste0(input,"cor_all.var_mgs.tsv"))
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   BMI = res[exposure=="BMI",])
  
  res.list = lapply(res.list,Clean.Correlation.Results)
  
  # Restricting results to mgs.m1
  res.list = lapply(res.list,function(x){
    x[x$MGS %in% mgs.m1,]
  })
  
  # From long to wide
  res.df = rbind(res.list[[1]],res.list[[2]],res.list[[3]])
  #a = c("cor.","p.value","q.value","N")
  #res.df = dcast(setDT(res.df), MGS~exposure,value.var=a)
  
  res.df <- res.df %>% pivot_wider(id_cols = MGS, names_from=exposure, 
                                   values_from= c("cor.","p.value","q.value","N"))
  
  # Save table 
  saveRDS(res.df,file="/home/baldanzi/Sleep_apnea/Results/table.res1.rds")
  

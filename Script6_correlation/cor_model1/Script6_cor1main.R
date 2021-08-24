# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-24

rm(list = ls())

# Loading packages 
pacman::p_load(data.table, ppcor, fastDummies, vegan)

# Output folders 
input1 = '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/'
output = "/home/baldanzi/Sleep_apnea/Results/"

# Importing data
pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

#Calculate Shannon diversity ####
a = grep("____",names(pheno),value=T) # vector with MGS names 
pheno[,shannon:=diversity(pheno[,a, with=F],index="shannon")]

# Removing rare MGS   
  # Calculating MGS prevalence 
  noms=grep("____",names(pheno),value=T)
# presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
  data_pa <- decostand(x = pheno[,noms,with=F], "pa")
# calculate sum per species
  data_sum <- data.frame(prevalence=apply(data_pa, 2, sum),
                       percentage = apply(data_pa,2,sum)/nrow(pheno))
  data_sum$MGS = rownames(data_sum)
#a = data_sum$MGS[data_sum$prevalence<5] #19 MGS are only present in less than 5 individuals
  a = data_sum$MGS[data_sum$percentage<.01] #383 MGS are present in less than 1% of individuals

  pheno <- pheno[ , -a, with=F] 
  
  
# Correlations 
    source(paste0(input1,"Spearman.correlation.function.R")) # Correlation function 
    source('cor_model1/Script6_cor1_AHI_MGS.R')
    source('cor_model1/Script6_cor1_T90_MGS.R')
    source('cor_model1/Script6_cor1_BMI_MGS.R')
  
  mgs.m1.filter001 <- unique(c(res.ahi$MGS, res.t90$MGS , res.bmi$MGS))
  saveRDS(mgs.m1.filter001, paste0(output,'mgs.m1.filter001.rds'))
  

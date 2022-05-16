# Script 6 - Sensitivity Analysis - 3

# MGS correlated with AHI and T90 but not to BMI in the fully adjusted model 
# are evaluated in a model adjusted for BMI 

rm(list = ls())
  
  pacman::p_load(data.table, ppcor, fastDummies,vegan, tidyverse)

  output = "/home/baldanzi/Sleep_apnea/Results/"

  input = "/home/baldanzi/Sleep_apnea/Results/"

# Sensitivity analysis 3
  # Add BMI to the full.model
  
  # Import MGS identified in model 2
  res.m2 <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  
  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  mgs.bmi <- res.m2[exposure =="BMI" & q.value<.05,MGS] # MGS correlated to BMI
  
  mgs.rel <- res.m2 %>% filter(exposure =="ahi" | exposure == "t90") %>% 
    filter(q.value<.05) %>% filter(!MGS %in% mgs.bmi) %>% select(MGS)
    
  mgs.rel <- unique(mgs.rel$MGS) # MGS identified in model 2 #

  # Importing data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  
  a= c("Sex", "metformin","hypermed","dyslipmed","ppi")
  pheno[,(a):=as.data.frame(data.matrix(data.frame(unclass(pheno[,a, with=F]))))]

#  Sens. Analysis 2
  #  adjust for model2 + sleeptime
  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  model2 <- c(model1,"metformin","hypermed","dyslipmed","ppi","Fibrer",
            "Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")
  SA3 <- c(model2,"BMI")

  # Prepare outcomes
  outcomes <-  mgs.rel

  #Spearman correlation function ####
  source("Spearman.correlation.function.R")

  # Correlation scripts   
  source("cor_SA3/Script6_corSA3_AHI_MGS.R")
  
  source("cor_SA3/Script6_corSA3_T90_MGS.R")
  
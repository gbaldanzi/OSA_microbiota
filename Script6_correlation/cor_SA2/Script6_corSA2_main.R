# Script 6 - Sensitivity Analysis - 2

rm(list = ls())
  
  pacman::p_load(data.table, ppcor, fastDummies,vegan)

  output = "/home/baldanzi/Sleep_apnea/Results/"

  input = "/home/baldanzi/Sleep_apnea/Results/"

# Sensitivity analysis 
# Remove users of ppi, metformin, hypertensive and cholesterol-lowering medication 
#  adjust for model1 + place birth + education + leisure PA + fiber + total energy intake 

  # Import MGS identified in model 1
  mgs.m1  = readRDS(paste0(input,'mgs.m1.rds'))

  # Importing data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  
  a= c("Sex", "metformin","hypermed","dyslipmed","ppi")
  pheno[,(a):=as.data.frame(data.matrix(data.frame(unclass(pheno[,a, with=F]))))]

#  Sens. Analysis 2
  #  adjust for model2 + sleeptime
  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  model2 <- c(model1,"metformin","hypermed","dyslipmed","ppi","Fibrer",
            "Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")
  SA2 <- c(model2,"sleeptime")

  # Prepare outcomes
  outcomes <-  mgs.m1

  #Spearman correlation function ####
  source("Spearman.correlation.function.R")

  # Correlation scripts   
  source("cor_SA2/Script6_corSA2_AHI_MGS.R")
  
  source("cor_SA2/Script6_corSA2_T90_MGS.R")
  
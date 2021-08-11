# Script 6 - Sensitivity Analysis - 2

rm(list = ls())

output = "/home/baldanzi/Sleep_apnea/Results/"

#  Sens. Analysis 2
#  adjust for model2 + sleeptime
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
model2 <- c(model1,"metformin","hypermed","dyslipmed","ppi","Fibrer",
            "Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")
SA2 <- c(model2,"sleeptime")


  # Loading packages 
  pacman::p_load(data.table, ppcor, fastDummies,vegan)

#Spearman correlation function ####
source("Spearman.correlation.function.R")

# Correlation scripts   
  source("cor_SA2/Script6_corSA2_AHI_MGS.R")
  
  source("cor_SA2/Script6_corSA2_T90_MGS.R")
  
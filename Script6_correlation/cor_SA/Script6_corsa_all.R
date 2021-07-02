# Script 6 - Sensitivity analysis 

rm(list = ls())

output = "/home/baldanzi/Sleep_apnea/Results/"

# Sensitivity analysis 
# Remove users of ppi, metformin, hypertensive and cholesterol-lowering medication 
#  adjust for model1 + place birth + education + leisure PA + fiber + total energy intake 

model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")

model2 <- c(model1,"Fibrer","Energi_kcal" ,"leisurePA", "educat","placebirth")

# Loading packages 
pacman::p_load(data.table, ppcor, fastDummies,vegan)

#Spearman correlation function ####
source("Spearman.correlation.function.R")

# Correlation scripts   
  source("cor_SA/Script6_corsa_AHI_MGS.R")
  
  source("cor_SA/Script6_corsa_BMI_MGS.R")
  
  source("cor_SA/Script6_corsa_T90_MGS.R")
  
# Script 6 - Model 2

rm(list = ls())

output = "/home/baldanzi/Sleep_apnea/Results/"

#Model 2
#  adjust for model1 + place birth + education + leisure PA + fiber + total energy intake +
#  diabetes + hypertension + dyslipidemia, medication 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
model2 <- c(model1,"metformin","hypermed","dyslipmed","ppi","Fibrer",
            "Energi_kcal" ,"leisurePA", "educat","placebirth")

# Loading packages 
pacman::p_load(data.table, ppcor, fastDummies,vegan)

#Spearman correlation function ####
source("Spearman.correlation.function.R")

# Correlation scripts   
  source("cor_model2/Script6_cor2_AHI_MGS.R")
  
  source("cor_model2/Script6_cor2_BMI_MGS.R")
  
  source("cor_model2/Script6_cor2_T90_MGS.R")
  
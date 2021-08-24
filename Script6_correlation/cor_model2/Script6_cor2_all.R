# Script 6 - Model 2

  rm(list = ls())

  output = "/home/baldanzi/Sleep_apnea/Results/"

  input = "/home/baldanzi/Sleep_apnea/Results/"

# Importing results 
  #res.ahi <- fread(paste0(input,"cor_ahi_mgs.tsv"))
  #res.bmi <- fread(paste0(input,"cor_BMI_mgs.tsv"))
  #res.t90 <- fread(paste0(input,"cor_t90_mgs.tsv"))
  #res.list = list(res.ahi,res.bmi,res.t90)
  
  res.ahi <- fread(paste0(input,"cor_ahi_mgs_filter001.tsv"))
  res.bmi <- fread(paste0(input,"cor_BMI_mgs_filter001.tsv"))
  res.t90 <- fread(paste0(input,"cor_t90_mgs_filter001.tsv"))

# filter MGS significant at the FDR p-value<0.05
  res.list <- lapply(res.list, function(x){x[q.value>=0.001, q.value:=round(q.value, digits = 3)]})
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})

  mgs.m1  = unique(unlist(mgs.fdr))

#Model 2
#  adjust for model1 + place birth + education + leisure PA + fiber + total energy intake +
#  diabetes + hypertension + dyslipidemia, medication 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
model2 <- c(model1,"metformin","hypermed","dyslipmed","ppi","Fibrer",
            "Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")

# Loading packages 
pacman::p_load(data.table, ppcor, fastDummies,vegan)

#Spearman correlation function ####
source("Spearman.correlation.function.R")

# Correlation scripts   
  source("cor_model2/Script6_cor2_AHI_MGS.R")
  
  source("cor_model2/Script6_cor2_BMI_MGS.R")
  
  source("cor_model2/Script6_cor2_T90_MGS.R")
  
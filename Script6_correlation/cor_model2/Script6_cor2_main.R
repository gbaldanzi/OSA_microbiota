# Script 6 - Model 2

  rm(list = ls())
  
  pacman::p_load(data.table, ppcor, fastDummies,vegan)

  output = "/home/baldanzi/Sleep_apnea/Results/"

  input = "/home/baldanzi/Sleep_apnea/Results/"
  
  
  # Import MGS identified in model 1
  mgs.m1  = readRDS(paste0(input,'mgs.m1.rds'))
  
  # Importing data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  
#Model 2
#  adjust for model1 + place birth + education + leisure PA + fiber + total energy intake +
#  diabetes + hypertension + dyslipidemia, medication 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
model2 <- c(model1,"metformin","hypermed","dyslipmed","ppi","Fibrer",
            "Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")


#Spearman correlation function ####
source("Spearman.correlation.function.R")

# Correlation scripts   
  source("cor_model2/Script6_cor2_AHI_MGS.R")

  source("cor_model2/Script6_cor2_T90_MGS.R")

  source("cor_model2/Script6_cor2_ODI_MGS.R")
  
  source("cor_model2/Script6_cor2_BMI_MGS.R")

  res <- rbind(res.ahi, res.t90, res.odi, res.bmi)
  
  # Merge taxonomic annotation 
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  setnames(taxonomy,"maintax_mgs","MGS")
  
  res <- merge(res, taxonomy, by="MGS", all.x=T)
  
  fwrite(res, file = paste0(output,"cor2_all.var_mgs.tsv"))
  
  res$q.value[res$q.value>=0.001] <- round(res$q.value[res$q.value>=0.001] , digits = 3)
  
  setDT(res)
  
  mgs.fdr.bmi = res[exposure=='BMI' & q.value<.05,MGS]
  mgs.fdr.ahi = res[exposure=="ahi" & q.value<.05 &!MGS %in% mgs.fdr.bmi,MGS]
  mgs.fdr.t90 = res[exposure=='t90' & q.value<.05 & !MGS %in% mgs.fdr.bmi,MGS]
  mgs.fdr.odi = res[exposure=='odi' & q.value<.05 & !MGS %in% mgs.fdr.bmi,MGS]
                     
  
  mgs.m2 <- list(mgs.fdr.ahi = mgs.fdr.ahi,
                 mgs.fdr.t90 = mgs.fdr.t90,
                 mgs.fdr.odi = mgs.fdr.odi) 
  
  
  saveRDS(mgs.m2, paste0(output,'mgs.m2.rds'))
 
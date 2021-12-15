# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Created: 2021-10-21

# This script was created to carry out the sensitivity analysis comparing the results 
# after exclusion of participants who have used antibiotic 

# In this script, we will only include the mgs that were FDR significant in the basic model 

  library(data.table)
  library(dplyr)
  library(tidyr)
  pacman::p_load(ppcor, fastDummies, vegan)
  
  
  input1 = '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/'
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output = "/home/baldanzi/Sleep_apnea/Results/"

  
  # Import data 
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

# Antibiotic data 
  atb = fread('/home/baldanzi/Datasets/Antibiotics/Processed/scapis_antibiotics_data_v1_20200804.tsv')
  atb <- atb[Site=="Site5",]
  atb_SCAPISid <- atb[Time_J01>-180,SUBJID]
  
  # Flagging individuals that used atb in the last 3 months. 
  pheno[!SCAPISid %in% atb_SCAPISid, atb6m:='no']  
  pheno[SCAPISid %in% atb_SCAPISid, atb6m:='yes']
  
  pheno <- pheno[atb6m=="no",]
  
  # Import results 
  
  #res <- fread(paste0(input,"cor_sa_atb3m_all.var_mgs.tsv"))
  
  #res[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  mgs.fdr  = readRDS(paste0(input,'mgs.m2.rds'))  # Signature MGSs 
  
  
  #Model 2
  #  adjust for model1 + place birth + education + leisure PA + fiber + total energy intake +
  #  diabetes + hypertension + dyslipidemia, medication 
  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  model2 <- c(model1,"metformin","hypermed","dyslipmed","ppi","Fibrer",
              "Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")
  
  # Correlations 
  source(paste0(input1,"Spearman.correlation.function.R")) # Correlation function 
  
  mgs.m1 = mgs.fdr$mgs.fdr.ahi
  
  source('cor_model2/Script6_cor2_AHI_MGS.R')
  
  mgs.m1 = mgs.fdr$mgs.fdr.t90
  
  source('cor_model2/Script6_cor2_T90_MGS.R')
  
  res <- rbind(res.ahi, res.t90)
  
  res$model= "sa_atb6m"
  
  names(res) = c("MGS", "exposure", "rho", "p.value", 
                 "N", "method", "covariates","q.value","model")
  
  fwrite(res, file = paste0(output,"cor_sa_atb6m_step2_all.var_mgs.tsv"))

  # Final table and venn diagram 
  message("Table of Results")
  source('cor_SA_atb/Script6_sa_atb_table.res2.R')
  message("Producing the Venn diagram")
  source('cor_SA_atb/Script6_sa_atb_venn2.R')

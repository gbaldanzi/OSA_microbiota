# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Created: 2021-10-21

# This script was created to carry out the sensitiviy analysis comparing the results 
# after exclusion of participants who have used antibiotic 

  library(data.table)
  library(dplyr)
  library(tidyr)
  pacman::p_load(ppcor, fastDummies, vegan)
  
  
  input1 = '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/'
  output = "/home/baldanzi/Sleep_apnea/Results/"

  
  # Import data 
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

# Antibiotic data 
  atb = fread('/home/baldanzi/Datasets/Antibiotics/Processed/scapis_antibiotics_data_v1_20200804.tsv')
  atb <- atb[Site=="Site5",]
  atb_SCAPISid <- atb[Time_J01>-90,SUBJID]
  
  # Flagging individuals that used atb in the last 3 months. 
  pheno[!SCAPISid %in% atb_SCAPISid, atb3m:='no']  
  pheno[SCAPISid %in% atb_SCAPISid, atb3m:='yes']
  
  # Remove rare MGS 
  mgs=grep("____",names(pheno[atb3m=='no',]),value=T)
  data_pa <- decostand(x = pheno[atb3m=='no',mgs,with=F], "pa")
  data_sum <- data.frame(prevalence=apply(data_pa, 2, sum))
  data_sum$MGS <- rownames(data_sum) 
  rarespecies <-  data_sum[data_sum$prevalence<5,"MGS"]
  
  pheno <- pheno[atb3m=='no' , -rarespecies, with=F] 
  
  # Correlations 
  source(paste0(input1,"Spearman.correlation.function.R")) # Correlation function 
  source('cor_model1/Script6_cor1_AHI_MGS.R')
  source('cor_model1/Script6_cor1_T90_MGS.R')
  source('cor_model1/Script6_cor1_BMI_MGS.R')
  
  res <- rbind(res.ahi, res.t90, res.bmi)
  
  res$model= "sa_atb3m"
  
  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                 "N", "method", "covariates","q.value","model")
  
  # Merging results with taxonomy information #### 
  
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  setnames(taxonomy,"maintax_mgs","MGS")
  
  res <- merge(res, taxonomy, by="MGS", all.x=T)
  
  fwrite(res, file = paste0(output,"cor_sa_atb3m_all.var_mgs.tsv"))



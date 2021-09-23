# Script 6 - Sensitivity analysis 

# Gabriel Baldanzi 
# 2021-09-02

 message("Sensitivity Analysis 1")

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
  
  # Remove medication users
  pheno <-  pheno[ppi == "no",] #60
  pheno <-  pheno[metformin == "no",] #57 
  pheno <-  pheno[hypermed == "no",] #593
  pheno <-  pheno[dyslipmed == "no",] # 239
  
  a= c("Sex")
  pheno[,(a):=as.data.frame(data.matrix(data.frame(unclass(pheno[,a, with=F]))))]

  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  model2 <- c(model1,"Fibrer","Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")
  
  outcomes <- mgs.m1


#Spearman correlation function ####
  source("/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/Spearman.correlation.function.R")

# Correlation scripts   
  source("cor_SA/Script6_corsa_AHI_MGS.R")
  
  source("cor_SA/Script6_corsa_T90_MGS.R")
  
  source("cor_SA/Script6_corsa_BMI_MGS.R")
  
  
  res <- rbind(res.ahi, res.t90, res.bmi)
  
  # Merge taxonomic annotation 
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  setnames(taxonomy,"maintax_mgs","MGS")
  
  res <- merge(res, taxonomy, by="MGS", all.x=T)
  
  fwrite(res, file = paste0(output,"corsa_all.var_mgs.tsv"))
  

# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Version 1: 2021-09-17
# Last Update: 2021-09-17

  rm(list = ls())

  pacman::p_load(data.table,ppcor, fastDummies, vegan)

# Defining output folders 
  output = "/home/baldanzi/Sleep_apnea/Results/"

# Import data
  pheno.self <- readRDS('/home/baldanzi/Datasets/sleep_SCAPIS/self.apnea.rds')

#Spearman correlation function ####
source('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/Spearman.correlation.function.R')

#-----------------------------------------------------------------------------#
# Models ####

  #Covariates 
  # model 1 : adjust for age + sex + alcohol + smoking + plate + received 
  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate", "shannon")
  # model 2 = model 1 + BMI 
  model2 <-  c(model1,"BMI")
  # model 3 = model 2 + fiber intake + Energy intake + physical activity + education + country of birth + ppi + metformin +  antihypertensive + cholesterol-lowering 
  #model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", "visit.month", "metformin","hypermed","dyslipmed","ppi")

  #listmodels=list(model1,model2,model3,SA2)
  

#-----------------------------------------------------------------------------#
# Correlation between AHI and Shannon ####
  
  # Prepare data 
  dades = copy(pheno.self[!is.na(self.apnea)])
  
  dades[self.apnea=='YES',self.apnea:=1]
  dades[self.apnea=='NO',self.apnea:=0]
  dades[,self.apnea:=as.numeric(self.apnea)]
  
  # Transforming two-level factor variables into numeric variables 
  a= c("Sex")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
  
# Preparing exposure and  outcomes
  exposure="self.apnea"
  
  # Outcome
  res.m2 <-  fread("/home/baldanzi/Sleep_apnea/Results/cor2_all.var_mgs.tsv")
  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  mgs.bmi <- res.m2[q.value<0.05 & exposure=="BMI",MGS]
  mgs.fdr <- unique(res.m2[q.value<0.05 & exposure %in% c("ahi","t90"),MGS])
  mgs.fdr <- mgs.fdr[!mgs.fdr %in% mgs.bmi] # The identified MGS #
  
  outcomes <- mgs.fdr

# Run Spearman correlation for the models.
  res <-  spearman.function(x1=outcomes,x2=exposure,covari = model1,data = dades)

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                       "N", "method", "covariates","q.value")

  fwrite(res, file=paste0(output,"cor_self.apnea_mgs.tsv"))
  

# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: March, 2021
# Update: July 8, 2021

# This code will investigate the alpha diversity in relation AHI, OSAcat, and T90%.  
# Four models  and a sensitivity analysis will be used 
# Results will be saved in three files depending on the exposure (AHI,OSAcat,T90%)

  rm(list = ls())

  pacman::p_load(data.table,ppcor, fastDummies, vegan)

# Defining output folders 
  output = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
  valid.ahi <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/validsleep_MGS.shannon_Upp.rds") 
  
  # Transforming into factor variable
  valid.ahi[,plate:=as.factor(valid.ahi$plate)]


#Spearman correlation function ####
source('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/Spearman.correlation.function.R')

#-----------------------------------------------------------------------------#
# Models ####

# model 1 : adjust for age + sex + alcohol + smoking + plate
  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate")
# model 2 = model 1 + BMI 
  model2 <-  c(model1,"BMI")
# model 3 = model 2 + fiber intake +energy intake+ physical activity + education + country of birth + visit month + medication
  model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth","visit.month", 
               "metformin","hypermed","dyslipmed","ppi")
# SA = model3 but removed medication users 
  SA <- c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth","visit.month")
  
  # SA2 = model 3 but added sleep duration 
  SA2 <- c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth","visit.month","sleeptime")
  
  listmodels=list(model1,model2,model3,SA2)
  
  
# OLD model 3 = model 2 + diabetes + hypertension + dyslipidemia, medication 
# OLD model3 = c(model2,"diabd","hypertension","dyslipidemia","diabmed","hypermed","dyslipmed","ppi")

#-----------------------------------------------------------------------------#
# Correlation between AHI and Shannon ####
  
  # Prepare data 
  dades = copy(valid.ahi)
  
  # Transforming two-level factor variables into numeric variables 
  a= c("Sex","hypermed","dyslipmed","ppi", "metformin")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
  
# Preparing exposure and  outcomes
  exposure="ahi"
  outcomes="shannon"

# Run Spearman correlation for the models.
  res.alpha = t(sapply(listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)}))

  res.alpha <-  as.data.frame(res.alpha)
  res.alpha$model <-  c("model1", "model2", "model3","SA2")


#----------------------------------------------------------------------------#
# Sensitivity analysis 1
  
  # Prepare data 
  dades.sa = copy(valid.ahi)
  dades.sa <- dades.sa[which(ppi=="no" & metformin=="no" & hypermed == "no" & dyslipmed =="no"), ]
  
  # Transforming two-level factor variables into numeric variables 
  a= c("Sex")
  dades.sa[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades.sa[,a, with=F]))))]

  # Spearman correlation 
  res.sa = spearman.function(x1=outcomes,x2=exposure,covari = SA,data = dades.sa)
  res.sa$model= "SA"

  # Naming columns
  res.alpha = rbind(res.alpha,res.sa)
  
  names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","model")

  fwrite(res.alpha, file = paste0(output,"cor_ahi_alpha.tsv"), sep="\t")

#-----------------------------------------------------------------------------#

# Correlation between AHI and Shannon by BMI group ####

# the correlation between AHI and shannon will be separately by the 3 BMI groups:
# <25, >=25 & <30, and >=30.
  
  setDT(dades)
  setDT(dades.sa)

  # By BMI group 
  for(group in unique(valid.ahi[,BMIcat])){
    dades2 = dades[BMIcat==group,]


  # Run Spearman correlation for the models.
  res.alpha = sapply(listmodels,function(mod){
     spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades2)})
  res.alpha <- as.data.frame(t(res.alpha))
  
  res.alpha$model = c("model1", "model2", "model3","SA2")
  
  res.alpha$bmi= group


  
#----------------------------------------------------------------------------#
  print("Sensitivity analysis")
  # Sensitivity analysis 

  # BMI group
  dades2 <- dades.sa[dades.sa$BMIcat==group,]

  # Spearman correlation 
  res.sa = spearman.function(x1=outcomes,x2=exposure,covari = SA,data = dades2)

  res.sa$model= "SA"
  res.sa$bmi= group

  res.alpha = rbind(res.alpha,res.sa)
  
  #naming coluns
  names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","model","bmi")

  fwrite(res.alpha, file = paste0(output,"cor_ahi_alpha_bmi",group,".tsv"), sep="\t")
  print(res.alpha)
  }


# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: March, 2021
# Update: July 8, 2021

# This code will investigate the alpha diversity in relation AHI, OSAcat, and T90%.  
# Four models  and a sensitivity analysis will be used 
# Results will be saved in three files depending on the exposure (AHI,OSAcat,T90%)


# Loading packages 
pacman::p_load(data.table,ppcor, fastDummies, vegan, parallel)

# Cleaning the environment 
rm(list = ls())

# Defining output folders 
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.ahi <- readRDS("validsleep_MGS.shannon_Upp.rds") 
setnames(valid.ahi, "pob", "placebirth")

# Transforming two-level factor variables into numeric variables 
a= c("Sex","hypermed","dyslipmed","ppi", "metformin")
valid.ahi[,(a):=as.data.frame(data.matrix(data.frame(unclass(valid.ahi[,a, with=F]))))]

# Transforming factor variables 
valid.ahi[,plate:=as.factor(valid.ahi$plate)]

# Prepare data 
  dades = copy(valid.ahi)
  
  # For Sensitivity analysis - remove individuals who use medication 
  dades.sa = copy(valid.ahi)
  dades.sa <-  dades.sa[dades.sa$ppi == 1,] #60
  dades.sa <-  dades[dades.sa$metformin == 1,] #57 
  dades.sa <-  dades[dades.sa$hypermed == 1,] #593
  dades.sa <-  dades[dades.sa$dyslipmed == 1,] # 239
  nrow(dades.sa) #2318

#Spearman correlation function ####
source(Spearman.correlation.function.R)

#-----------------------------------------------------------------------------#
# Models ####

# model 1 : adjust for age + sex + alcohol + smoking + plate
  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate")
# model 2 = model 1 + BMI 
  model2 <-  c(model1,"BMI")
# model 3 = model 2 + fiber intake +energy intake+ physical activity + education + country of birth + PPI + metformin + anti-hypertensive + cholesterol-lowering
  model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", 
               "metformin","hypermed","dyslipmed","ppi")
# SA = model3 but removed medication users 
  SA <- c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth")
  
  
  listmodels=list(model1,model2,model3)
  
  
# OLD model 3 = model 2 + diabetes + hypertension + dyslipidemia, medication 
# OLD model3 = c(model2,"diabd","hypertension","dyslipidemia","diabmed","hypermed","dyslipmed","ppi")

#-----------------------------------------------------------------------------#
# Correlation between AHI and Shannon ####

# Preparing exposure and  outcomes
  exposure="ahi"
  outcomes="shannon"

# Run Spearman correlation for the models.
res.alpha = lapply(listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)})
names(res.alpha) = c("model1", "model2", "model3")

res.model1 = res.alpha[[1]]
res.model2 = res.alpha[[2]]
res.model3 = res.alpha[[3]]

res.model1$model= "model1"
res.model2$model= "model2"
res.model3$model= "model3"

res.alpha = as.data.table(rbind(res.model1, res.model2, res.model3, res.model4))

names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","model")
#----------------------------------------------------------------------------#
# Sensitivity analysis 

# Spearman correlation 
res.sa = spearman.function(x1=outcomes,x2=exposure,covari = SA,data = dades.sa)

res.sa$model= "SA"

#naming coluns
names(res.sa) <- c("MGS", "exposure", "cor.coefficient", "p.value", 
                "N", "method", "covariates","model")

res.alpha = rbind(res.alpha,res.sa)

fwrite(res.alpha, file = paste0(output,"cor_ahi_alpha.tsv"), sep="\t")

#-----------------------------------------------------------------------------#

# Correlation between AHI and Shannon by BMI group ####

# the correlation between AHI and shannon will be separately by the 3 BMI groups:
# <25, >=25 & <30, and >=30.

# By BMI group 
for(group in unique(valid.ahi[,BMIcat])){
  dades2 = dades[BMIcat==group,]


# Run Spearman correlation for the models.
res.alpha = lapply(listmodels,function(mod){
     spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades2)})
names(res.alpha) = c("model1", "model2", "model3")
  
res.model1 = res.alpha[[1]]
res.model2 = res.alpha[[2]]
res.model3 = res.alpha[[3]]


res.model1$model= "model1"
res.model2$model= "model2"
res.model3$model= "model3"

res.model1$bmi= group
res.model2$bmi= group
res.model3$bmi= group

res.alpha = as.data.table(rbind(res.model1, res.model2, res.model3))

names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","model","bmi")
print(res.alpha)

#----------------------------------------------------------------------------#
print("Sensitivity analysis")
# Sensitivity analysis 

# BMI group
dades2 <- dades.sa[dades.sa$BMIcat==group,]

# Spearman correlation 
res.sa = spearman.function(x1=outcomes,x2=exposure,covari = SA,data = dades2)

res.sa$model= "SA"
res.sa$bmi= group

#naming coluns
names(res.sa) <- c("MGS", "exposure", "cor.coefficient", "p.value", 
                   "N", "method", "covariates","model","bmi")

res.alpha = rbind(res.alpha,res.sa)

fwrite(res.alpha, file = paste0(output,"cor_ahi_alpha_bmi",group,".tsv"), sep="\t")
}
#-----------------------------------------------------------------------------#

# Correlation between T90% and Shannon ####

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.t90 <- readRDS("valid.t90_MGS.shannon_Upp.rds")
setnames(valid.t90, "pob", "placebirth")


# Transforming two-level factor variables into numeric variables 
a= c("Sex","hypermed","dyslipmed","ppi", "metformin")
valid.t90[,(a):=as.data.frame(data.matrix(data.frame(unclass(valid.t90[,a, with=F]))))]

# Transforming factor variables 
valid.t90[,plate:=as.factor(valid.t90$plate)]

# Preparing exposure, outcomes, and covariates
exposure="t90"
outcomes="shannon"

# Run Spearman correlation for the models.
res.alpha = lapply(listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)})
names(res.alpha) = c("model1", "model2", "model3")

res.model1 = res.alpha[[1]]
res.model2 = res.alpha[[2]]
res.model3 = res.alpha[[3]]

res.model1$model= "model1"
res.model2$model= "model2"
res.model3$model= "model3"

res.alpha = as.data.table(rbind(res.model1, res.model2, res.model3,res.model4))

names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","model")
#----------------------------------------------------------------------------#
# Sensitivity analysis - remove individuals who use medication 
dades.sa = copy(valid.t90)
dades.sa <-  dades.sa[dades.sa$ppi == 1,] #60
dades.sa <-  dades.sa[dades.sa$metformin ==1,] #57 
dades.sa <-  dades.sa[dades.sa$hypermed == 1,] #593
dades.sa <-  dades.sa[dades.sa$dyslipmed == 1,] # 239
nrow(dades) #2602

# Spearman correlation 
res.sa = spearman.function(x1=outcomes,x2=exposure,covari = model3,data = dades.sa)

res.sa$model= "SA"

#naming coluns
names(res.sa) <- c("MGS", "exposure", "cor.coefficient", "p.value", 
                   "N", "method", "covariates","model")

res.alpha = rbind(res.alpha,res.sa)

#Saving results 
fwrite(res.alpha, file = paste0(output,"cor_t90_alpha.tsv"), sep="\t")

#-----------------------------------------------------------------------------#

# Correlation between T90% and Shannon by BMI group ####

# the correlation between AHI and shannon will be separetely by the 3 BMI groups:
# <25, >=25 & <30, and >=30.


#By BMI group
for(group in unique(valid.ahi[,BMIcat])){

dades2 = dades[BMIcat==group,]
  
# Run Spearman correlation for the models.
res.alpha = lapply(listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades2)})
names(res.alpha) = c("model1", "model2", "model3")

res.model1 = res.alpha[[1]]
res.model2 = res.alpha[[2]]
res.model3 = res.alpha[[3]]

res.model1$model= "model1"
res.model2$model= "model2"
res.model3$model= "model3"

res.model1$bmi= group
res.model2$bmi= group
res.model3$bmi= group

res.alpha = as.data.table(rbind(res.model1, res.model2, res.model3))

names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","model","bmi")
#----------------------------------------------------------------------------#
print("Sensitivity analysis (no medication users) correlation T90% and Shannon Index")
# Sensitivity analysis - remove individuals who use medication 

dades2 = dades.sa[BMIcat==group,]

# Spearman correlation 
res.sa = spearman.function(x1=outcomes,x2=exposure,covari = model3,data = dades2)

res.sa$model= "SA"
res.sa$bmi= group

#naming coluns
names(res.sa) <- c("MGS", "exposure", "cor.coefficient", "p.value", 
                   "N", "method", "covariates","model","bmi")

res.alpha = rbind(res.alpha,res.sa)

#Saving results 
fwrite(res.alpha, file = paste0(output,"cor_t90_alpha_bmi",group,".tsv"), sep="\t")
}



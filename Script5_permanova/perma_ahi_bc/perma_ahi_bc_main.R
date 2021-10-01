# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: April 2021
# Update: - Jul 8, 2021: merged models 3 and 4, including medication with other covariates. 

# This code will investigate the beta-diversity (Bray Curtis Dissimilarity and 
# Aitchison distance) in relation to AHI and T90% in 3 different models. 
# Analysis are run using a PERMANOVA approach

# Results saved at the folder: "/home/baldanzi/Sleep_apnea/Results/"
# File model 1 - permanova_model1.tsv
# File model 2 - permanova_model2.tsv
# File model 3 - permanova_model3.tsv
# File SA - permanova_SA.tsv

# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)

rm(list = ls())
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
valid.ahi <- pheno[valid.ahi=='yes',]

# Importing BC matrix 
BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv', header=T, sep = ',')
BC = as.matrix(BC)
row.names(BC) = colnames(BC)

source('permanova.fun.R')

# Transforming two-level factor variables into numeric variables 
dades = copy(valid.ahi)
a= c("Sex","ppi","metformin","hypermed","dyslipmed")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]


# Making sure that BC and dades have the same order of observations 
dades = dades[match(rownames(BC),dades$SCAPISid),]

# Outcome - character name (length=1) with matrix distance 
outc = "BC"

# Main Exposure - character name (length=1)
expo = "ahi"

#Covariates 
  # model 1 : adjust for age + sex + alcohol + smoking + plate + received 
  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate")
  # model 2 = model 1 + BMI 
  model2 <-  c(model1,"BMI")
  # model 3 = model 2 + fiber intake + Energy intake + physical activity + education + country of birth + ppi + metformin +  antihypertensive + cholesterol-lowering 
  model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", "visit.month", "metformin","hypermed","dyslipmed","ppi")
  # SA = remove medication users 
  SA <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", "visit.month")
  # SA2 = model 2 + fiber intake + Energy intake + physical activity + education + country of birth + ppi + metformin +  antihypertensive + cholesterol-lowering 
  SA2 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", "visit.month", "sleeptime", "metformin","hypermed","dyslipmed","ppi")
  
  
  
# Runing PERMANOVA in parallel ####
  source('perma_ahi_bc/perma_model1.R')
  
  source('perma_ahi_bc/perma_model2.R')

  source('perma_ahi_bc/perma_model3.R')

  source('perma_ahi_bc/perma_SA.R')
  
  source('perma_ahi_bc/perma_SA2.R')

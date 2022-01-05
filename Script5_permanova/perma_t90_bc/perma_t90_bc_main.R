# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: April 2021
# Update: 2021-12-09 

# This code will investigate the beta-diversity (Bray Curtis Dissimilarity)
# in relation to T90% severity categories 
# Analysis are run using a PERMANOVA approach

#parallel nodes (1 for each model). Therefore, this code requires at least 3 nodes. 

# Results saved at the folder: "/home/baldanzi/Sleep_apnea/Results/"
# File model 3 - permanova_model3_t90_bc.tsv


# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)

rm(list = ls())
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing data

  pheno <- readRDS(file="/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")


  # Importing BC matrix 
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/T90.BCmatrix.csv', header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)

  source('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/permanova.fun.R')

  # Transforming two-level factor variables into numeric variables 
  dades = copy(pheno[ valid.t90 == 'yes', ])
  a= c("Sex","ppi","metformin","hypermed","dyslipmed")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]


  # Making sure that BC and dades have the same order of observations 
  dades = dades[match(rownames(BC),dades$SCAPISid),]

  # Outcome - character name (length=1) with matrix distance 
  outc = "BC"

  # Main Exposure - character name (length=1)
  expo = "t90cat"

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
  
  

# Runing PERMANOVA ####
  
  #message('perma_t90_bc/perma_model1.R')
  #source('perma_t90_bc/perma_model1.R')

  #message('perma_t90_bc/perma_model2.R')
  #source('perma_t90_bc/perma_model2.R')

  message('perma_t90_bc/perma_model3.R')
  source('perma_t90_bc/perma_model3.R')
  
  #message('perma_t90_bc/perma_SA2.R')
  #source('perma_t90_bc/perma_SA2.R')

  #message('perma_t90_bc/perma_SA.R')
  #source('perma_t90_bc/perma_SA.R')
  



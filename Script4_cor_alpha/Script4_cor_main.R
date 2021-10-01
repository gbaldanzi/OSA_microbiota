# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlation between Alpha diversity and the phenotypes T90 and AHI 
# Update: Sep 30, 2021

  pacman::p_load(data.table,ppcor, fastDummies, vegan, ggplot2, grid)

# Defining output folders 
  output = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")


#Spearman correlation function ####
source('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/Spearman.correlation.function.R')

#-----------------------------------------------------------------------------#
# Models ####

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
  
  
  listmodels=list(model1,model2,model3,SA2)
  
# Run correlation analysis 
  
  message("Correlation AHI and Shannon index")
  source('Script4_cor_alpha/Script4_cor_ahi_alpha.R')

  message("Correlation T90 and Shannon index")
  source('Script4_cor_alpha/Script4_cor_t90_alpha.R')
  
  message("Correlation AHI and Shannon index - Plots")
  source('Script4_cor_alpha/Script4_alpha_ahi_plots.R')
  
  message("Correlation T90 and Shannon index - Plots")
  source('Script4_cor_alpha/Script4_alpha_t90_plots.R')
  
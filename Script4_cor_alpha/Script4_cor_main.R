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
  # basic.model : adjust for age + sex + alcohol + smoking + plate + received 
  basic.model <-   c("age", "Sex", "Alkohol","smokestatus","plate")
  # full model = basic model + fiber intake + Energy intake + physical activity + education + country of birth + ppi + metformin +  antihypertensive + cholesterol-lowering 
  full.model <-  c(basic.model, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", "visit.month", "metformin","hypermed","dyslipmed","ppi")
  # full.model_bmi = full model + BMI 
  full.model_bmi = c(full.model, "BMI")
  
  #listmodels=list(model1,model2,model3,SA2)
  listmodels = list(basic.model, full.model, full.model_bmi)
  names(listmodels) <- c("basic model", "full model", "full model + BMI")
  
# Run correlation analysis 
  
  message("Correlation AHI and Shannon index")
  source('Script4_cor_alpha/Script4_cor_ahi_alpha.R')
  
  message("Correlation ODI and Shannon index")
  source('Script4_cor_alpha/Script4_cor_t90_alpha.R')

  message("Correlation T90 and Shannon index")
  source('Script4_cor_alpha/Script4_cor_odi_alpha.R')
  
  res.alpha <- rbind(res.alpha.ahi, res.alpha.t90, res.alpha.odi)
  
  fwrite(res.alpha, file= paste0(output, "cor_all.var_alpha.tsv"))
  
  #message("Correlation AHI and Shannon index - Plots")
  #source('Script4_cor_alpha/Script4_alpha_ahi_plots.R')
  
  #message("Correlation T90 and Shannon index - Plots")
  #source('Script4_cor_alpha/Script4_alpha_t90_plots.R')
  
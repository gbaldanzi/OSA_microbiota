# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlation with species abundance

# This set of scripts will investigate the association between the relative 
# abundance of gut microbiota species and four phenotypes (AHI, T90, ODI, BMI)

  # Loading packages 
  library(data.table)
  library(tidyr)
  library(dplyr)

  # Folders 
  results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
  work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'
  
  # Import data
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))
  
  # Models
  main.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  main.model.BMI <- c(main.model, "BMI")
  extended.model <- c(main.model.BMI,"Fibrer","Energi_kcal" ,"leisurePA", 
                  "educat","placebirth","visit.month")
  medication.model <- c(main.model.BMI, "metformin","hypermed", "dyslipmed", "ppi")
  
  # Outcomes and exposures
  outcomes  <-  grep("HG3A",names(pheno),value=T)
  exposures <- c("ahi","t90","odi")
  
  # Functions 
  source("0_functions/Spearman.correlation.function.R") # Correlation function 
  #source("0_functions/Script6.Functions.R")


  message("Main Model without BMI")
  source("3_correlation/cor_main.model.R")
  
  message("Main Model with BMI")
  source("3_correlation/cor_main.model.BMI.R")
  
  message("Extended Model")
  source("3_correlation/cor_extended.model.R")
  
  message("Medication Model")
  source("3_correlation/cor_medication.model.R")
  
  message("Sensitivity Analysis - Atb")
  source("3_correlation/cor_sa_atb.R")
  
  message("Sensitivity Analysis - Lund Disease")
  source("3_correlation/cor_sa_lungdisease.R")

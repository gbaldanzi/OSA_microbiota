# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlation with species abundance

# This set of scripts will investigate the association between the OSA parameters (AHI, T90, ODI)
# and the relative abundance of gut microbiota species 

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
  main.model <-   c("age", "Sex", "Alkohol","smokestatus","plate")
  main.model.BMI <- c(main.model, "BMI")
  extended.model <- c(main.model.BMI,"Fibrer","Energi_kcal" ,"leisurePA", 
                  "educat","placebirth","visit.month")
  medication.model <- c(extended.model, "metformin","hypermed", "dyslipmed", "ppi")
  
  # Outcomes and exposures
  outcomes  <-  grep("HG3A",names(pheno),value=T)
  exposures <- c("ahi","t90","odi")
  
  # Function
  source("0_functions/Spearman.correlation.function.R") # Correlation function 

  # Main model excluded of BMI
  source("3_correlation/cor_main.model.R")
  
  # Complete main model 
  source("3_correlation/cor_main.model.BMI.R")
  
  # Extended model 
  source("3_correlation/cor_extended.model.R")
  
  # Species associated with AHI, T90, or ODI in the extended model are investigated 
  # in three sensitivity analyses
  
  ## Sensitiviy analysis 1 - medication model 
  source("3_correlation/cor_medication.model.R")
  
  ## Sensitivity analysis 2 - removed antibiotic users 
  source("3_correlation/cor_sa_atb.R")
  
  ## Sensitivity analysis 3 - removed participants with self-reported lung disease
  source("3_correlation/cor_sa_lungdisease.R")

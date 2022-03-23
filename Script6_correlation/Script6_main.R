# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script 6 Main - Correlation with species abundance

# This set of scripts will investigate the association between the relative 
# abundance of gut microbiota species and four phenotypes (AHI, T90, ODI, BMI)

  # Loading packages 
  library(data.table)
  library(tidyr)
  library(dplyr)

  # Folders 
  input = "/home/baldanzi/Datasets/sleep_SCAPIS/"
  results.folder = "/home/baldanzi/Sleep_apnea/Results/"

  # Importing data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

  # Models
  main.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  main.model.BMI <- c(main.model, "BMI")
  extended.model <- c(main.model.BMI,"Fibrer","Energi_kcal" ,"leisurePA", 
                  "educat","placebirth","visit.month")
  medication.model <- c(main.model.BMI, "metformin","hypermed", "dyslipmed", "ppi")
  
  # Outcomes and exposures
  outcomes  <-  grep("____",names(pheno),value=T)
  exposures <- c("ahi","t90","odi")
  
  # Functions 
  source("Script0_functions/Spearman.correlation.function.R") # Correlation function 
  source("Script0_functions/Script6.Functions.R")

# Import data 

  message("Main Model without BMI")
  source("Script6_correlation/cor_main.model.R")
  
  message("Main Model with BMI")
  source("Script6_correlation/cor_main.model.BMI.R")
  
  message("Extended Model")
  source("Script6_correlation/cor_extended.model.R")
  
  message("Medication Model")
  source("Script6_correlation/cor_medication.model.R")
  
  message("Sensitivity Analysis - Atb")
  source("Script6_correlation/cor_sa_atb.R")
  
  message("Creating Venn Diagram")
  source("Script6_correlation/Script6.venn.R")
  
  #message("Sensitivity Analysis 1")
  #source("cor_sa_med.R")
  

  
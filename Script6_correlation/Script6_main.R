# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script 6 Main - Correlation with species abundance

# This set of scripts will investigate the association between the relative 
# abundance of gut microbiota species and four phenotypes (AHI, T90, ODI, BMI)

  # Loading packages 
  library(data.table)
  library(vegan)
  library(tidyr)
  library(dplyr)

  # Folders 
  input = "/home/baldanzi/Datasets/sleep_SCAPIS/"
  results.folder = "/home/baldanzi/Sleep_apnea/Results/"

  # Importing data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

  # Models
  basic.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  full.model <- c(basic.model,"BMI","Fibrer","Energi_kcal" ,"leisurePA", 
                  "educat","placebirth","visit.month","metformin","hypermed",
                  "dyslipmed","ppi")
  
  # Functions 
  source("Script0_functions/Spearman.correlation.function.R") # Correlation function 
  source("Script0_functions/Script6.Functions.R")

# Import data 

  message("Basic Model")
  source("cor_basic.model.R")
  
  message("Full Model")
  source("cor_full.model.R")
  
  message("Creating Venn Diagram")
  source("Script6.venn.R")
  
  message("Sensitivity Analysis 1")
  source("cor_sa_med.R")
  
  message("Sensitivity Analysis 2")
  source("cor_sa_atb.R")
  
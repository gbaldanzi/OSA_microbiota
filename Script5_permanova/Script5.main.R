# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script 5 - Main 

# Script to run the permutation analysis between OSA/AHI/T90 and BC

# This script in meant to run PERMANOVA analysis in parallel in 16 cores

# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)

  # Folders
  results.folder = "/home/baldanzi/Sleep_apnea/Results/"
  input = "/home/baldanzi/Datasets/sleep_SCAPIS/"

  # Import data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))


#Covariates 
  # Models
  basic.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","BMI")
  full.model <- c(basic.model,"Fibrer","Energi_kcal" ,"leisurePA", 
                  "educat","placebirth","visit.month")
  #medication.model <- c(basic.model, "metformin","hypermed", "dyslipmed", "ppi")

  # Permanova Function 
  source('Script0_functions/permanova.fun.R')

  
  # Analysis by exposure variable 
  source('Script5_permanova/perma_ahi_bc/perma_ahi_bc.R')

  source('Script5_permanova/perma_odi_bc/perma_odi_bc.R')

  source('Script5_permanova/perma_t90_bc/perma_t90_bc.R')


# Script to run the pairwise PERMANOVA for groups of AHI/T90/ODI 

# Every PERMANOVA analysis takes approx. 12h to run in 16 cores. 
# Thus, this script takes about 4 x 3 x 12h = 144h to finish, if ran sequentially 

  # Loading packages
  library(data.table)
  library(vegan)
  library(parallel)

  # Folders
  input = "/home/baldanzi/Datasets/sleep_SCAPIS/"
  results.folder = "/home/baldanzi/Sleep_apnea/Results/"

  
  # Import data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

  
  # Functions 
  source('0_functions/perma.pairwise.fun.R')

  
  #Covariates 
  main.model<-   c("age", "Sex", "Alkohol","smokestatus","plate", "BMI")
  
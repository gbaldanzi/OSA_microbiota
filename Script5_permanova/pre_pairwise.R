# Script 5 - Main 

# Script to run the pairwise permutation analysis for overall 
# microbiota differences between group of OSA/AHI/T90 

# Every PERMANOVA analysis takes approx. 12h to run in 16 cores. 
# Thus, this script takes about 4 x 3 x 12h = 144h to finish. 

  # Loading packages 
  pacman::p_load(data.table, vegan ,parallel)

  # Folders
  output = "/home/baldanzi/Sleep_apnea/Results/"
  input = "/home/baldanzi/Datasets/sleep_SCAPIS/"

  
  # Import data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

  
  # Functions 
  source('Script0_functions/perma.pairwise.fun.R')

  
  #Covariates 
  basic.model<-   c("age", "Sex", "Alkohol","smokestatus","plate")
  full.model <-  c(basic.model,"BMI", "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", "visit.month", "metformin","hypermed","dyslipmed","ppi")
  
 

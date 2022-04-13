# Script to run the pairwise PERMANOVA for groups of AHI/T90/ODI 

# Every PERMANOVA analysis takes approx. 12h to run in 16 cores. 
# Thus, this script takes about 4 x 3 x 12h = 144h to finish, if ran sequentially 

  # Loading packages
  library(data.table)
  library(vegan)
  library(parallel)

  # Folders
  results.folder = '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'

  work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'
  
  # Import data
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))

  
  # Functions 
  source('0_functions/perma.pairwise.fun.R')

  #Covariates 
  main.model<-   c("age", "Sex", "Alkohol","smokestatus","plate", "BMI")
  
  # Importing BC matrix 
  BC = fread(paste0(work,'BCmatrix.tsv'))
  BC_rownames <- BC$rownames
  BC$rownames <- NULL
  
  BC <-  as.matrix(BC)
  colnames(BC) <- rownames(BC) <- BC_rownames
  

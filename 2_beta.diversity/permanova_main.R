
# Script to run the permutation analysis between AHI/T90/ODI and BC

# This script is meant to run PERMANOVA analysis in parallel in 16 cores

  # Loading packages 
  library(data.table)
  library(vegan)
  library(parallel)

  # Folders
  results.folder = '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
  work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'
  
  # Import data
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))
  

  # Models
  main.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","BMI")
  extended.model <- c(main.model,"Fibrer","Energi_kcal" ,"leisurePA", 
                  "educat","placebirth","visit.month")

  # Permanova Function 
  source('0_functions/permanova.fun.R')
  
  # Importing BC matrix 
  BC = fread(paste0(work,'BCmatrix.tsv'))
  BC_rownames <- BC$rownames
  BC$rownames <- NULL
  
  BC <-  as.matrix(BC)
  colnames(BC) <- rownames(BC) <- BC_rownames

  
  # Analysis by exposure variable 
  
  # By AHI severity groups ####
  expo = "OSAcat"
  pheno.ahi <- pheno[valid.ahi=="yes",]
  
  # Run PERMANOVA in parallel 
  
  # Main model 
  res1 <- Permanova.parallel.FUN(y = BC, exposure=expo, 
                                data = pheno.ahi , 
                                model = main.model, nod=16)
  
  fwrite(res1, file = paste0(results.folder,"permanova_main.model_ahi.tsv"))
  
  # Extended model
  res2 <- Permanova.parallel.FUN(y = BC, exposure=expo, 
                                data = pheno.ahi, 
                                model = extended.model, nod=16)
  
  fwrite(res2, file = paste0(results.folder,"permanova_extended.model_ahi.tsv"))
  
  
  
  # By T90 severity groups ####
  expo = "t90cat"
  pheno.t90 <- pheno[valid.t90=="yes",]
  
  # Run PERMANOVA in parallel 
  
  # Main model 
  res3 <- Permanova.parallel.FUN(y = BC, exposure=expo, 
                                data = pheno.t90, 
                                model = main.model, nod=16)
  
  fwrite(res3, file = paste0(results.folder,"permanova_main.model_t90.tsv"))
  
  # Extended model
  res4 <- Permanova.parallel.FUN(y = BC, exposure=expo, 
                                data = pheno.t90, 
                                model = extended.model, nod=16)
  
  fwrite(res4, file = paste0(results.folder,"permanova_extended.model_t90.tsv"))
  
  
  
  # By ODI severity groups ####
  expo = "odicat"
  
  # Runing PERMANOVA in parallel 
  
  # Main model 
  res5 <- Permanova.parallel.FUN(y = BC, exposure=expo, 
                                data = pheno.t90, 
                                model = main.model, nod=16)
  
  fwrite(res5, file = paste0(results.folder,"permanova_main.model_odi.tsv"))
  
  # Extended model
  res6 <- Permanova.parallel.FUN(y = BC, exposure=expo, 
                                data = pheno.t90, 
                                model = extended.model, nod=16)
  
  fwrite(res6, file = paste0(results.folder,"permanova_extended.model_odi.tsv"))



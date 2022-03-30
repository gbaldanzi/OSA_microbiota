# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script will run the correlation between the GMM enriched among the T90 associations
# and the plasma metabolites. Next, it will run an enrichment analysis in the metabolites

# Last update: 2022-03-17

  # load packages
  library(data.table)

  # Load correlation function 
  source("Script0_functions/Spearman.correlation.function.R")

  results.folder <-  "/home/baldanzi/Sleep_apnea/Results/"

  # Import data set with GMM abundances 
  input <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work'
  phenofull <- readRDS(paste0(input,"phenotype_upp_malm.rds"))
  osa.gmm <- grep("MF0",names(phenofull),value=T)
  
  
  # Covariates 
  covariates <- c("age", "Sex", "shannon", "plate", "Site", "placebirth")
  
  
  # Correlation 
  exposures <- osa.gmm
  outcomes <- grep("MET_",names(phenofull),value = T)
  
  res <- lapply(exposures,spearman.function, 
                           x1=outcomes,
                           covari = covariates,
                data = phenofull[metabolon_data==T,])
  
  res <- do.call(rbind,res)
  
  setDT(res)
  setnames(res, c("x","exposure"), c("metabolites", "modules"))
  
  
  # Saving results 
  fwrite(res, file = paste0(results.folder,"cor_gmm_metabolites.tsv"))
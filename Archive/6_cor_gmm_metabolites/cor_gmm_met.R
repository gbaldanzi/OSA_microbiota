# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script will run the correlation between the GMM enriched among the T90 associations
# and the plasma metabolites. Next, it will run an enrichment analysis in the metabolites

  library(data.table)

  # Load correlation function 
  source("0_functions/Spearman.correlation.function.R")

  # Folders
  results.folder <-  "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/"
  input <- "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/data_processed/"

  # Import data set with GMM abundances for SCAPIS-Uppsala and SCAPIS-Malmo participants 
  phenofull <- readRDS(paste0(input,"phenotype_upp_malm.rds"))
  
  # GMMs (Pathways enriched among the positive T90-species correlations)
  osa.gmm <- grep("MF0",names(phenofull),value=T)
  
  
  # Covariates 
  covariates <- c("age", "Sex", "shannon", "site_plate", "placebirth")
  
  
  # Correlation 
  exposures <- osa.gmm # GMMs abundance
  outcomes <- grep("MET_",names(phenofull),value = T) # Metabolites 
  
  res <- lapply(exposures,
                FUN = spearman.function, 
                x1=outcomes,
                covari = covariates,
                data = phenofull[metabolon_data==T,])
  
  res <- do.call(rbind,res)
  
  setDT(res)
  setnames(res, c("x","exposure"), c("metabolites", "modules"))
  
  
  # Saving results 
  fwrite(res, file = paste0(results.folder,"cor_gmm_metabolites.tsv"))
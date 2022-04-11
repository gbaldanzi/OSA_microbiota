# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-22

# Correlations with gut microbiota species using the extended model


library(data.table)
library(tidyr)
library(dplyr)

# Folders 
input = "/home/baldanzi/Datasets/sleep_SCAPIS/"
results.folder = "/home/baldanzi/Sleep_apnea/Results/"

# Import MGS identified in model 1
mgs.m1  = readRDS(paste0(results.folder,'mgs.m1.rds'))

# Importing data
pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

# Models
sa.wo.shannon <- c("age", "Sex", "Alkohol","smokestatus","plate","BMI")

# Outcomes and exposures
outcomes  <-  mgs.m1
exposures <- c("t90","odi")
  

# Correlations

  res.sa.wo.shannon <- lapply(exposures,spearman.function, 
                          x1=outcomes,
                          covari = sa.wo.shannon,
                          data = pheno)

  res.sa.wo.shannon  <- do.call(rbind,res.sa.wo.shannon)
  
  setDT(res.sa.wo.shannon)
  setnames(res.sa.wo.shannon,"x","MGS")
  
  # Import main analysis results 
  res.main <- fread(paste0(results.folder,"cor.bmi_all.var_mgs.tsv"))
  
  res.sa.wo.shannon <- merge(res.sa.wo.shannon,res.main, 
                             by = c("MGS","exposure"),
                             suffixes = c("_sa","_main"),
                             all.x=T, all.y=F)
  
  
  # Save results for the extended model
  fwrite(res.sa.wo.shannon, file = paste0(results.folder,"cor_sa.wo.shannon.tsv"))
  
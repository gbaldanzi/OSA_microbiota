# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 
  library(vegan)
  library(data.table)

# This script will calculate alpha diversity

# Import data 
  pheno <- readRDS(file="/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

#Calculate Shannon diversity ####

  mgs.names = grep("____",names(pheno),value=T) # vector with MGS names 
  pheno[,shannon:=diversity(pheno[,mgs.names, with=F],index="shannon")]
  
# Save data

  saveRDS(pheno, file="/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

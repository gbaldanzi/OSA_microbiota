# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will calculate alpha diversity


# Import data 
  pheno <- readRDS(file="/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

#Calculate Shannon diversity ####

  a = grep("____",names(pheno),value=T) # vector with MGS names 
  pheno[,shannon:=diversity(pheno[,a, with=F],index="shannon")]
  
# Save data

  saveRDS(pheno, file="/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

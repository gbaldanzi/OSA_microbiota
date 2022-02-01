# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script will calculate alpha diversity (Shannon diversity index) for all participants

  library(vegan)
  library(data.table)

  # input = folder containing the dataset 
  input="/home/baldanzi/Datasets/sleep_SCAPIS/"
  
  # Import data 
  pheno <- readRDS(file=paste0(input,"pheno.MGS.Upp.rds"))

#Calculate Shannon diversity ####

  mgs.names = grep("____",names(pheno),value=T) # vector with MGS names 
  pheno[,shannon:=diversity(pheno[,mgs.names, with=F],index="shannon")]
  
# Save data

  saveRDS(pheno, file=paste0(input,"pheno.MGS.Upp.rds"))
  
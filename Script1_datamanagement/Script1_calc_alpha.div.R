# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script will calculate alpha diversity (Shannon diversity index) for all participants

  library(vegan)
  library(data.table)

  # input = folder containing the dataset 
  input="/home/baldanzi/Datasets/sleep_SCAPIS/"
  
  # Import data 
  pheno <- readRDS(file=paste0(input,"pheno.MGS.Upp.rds"))
  
# Removing the MGS that are present in less than 1% of the individuals. 
  
  # Calculating MGS prevalence 
  species.names=grep("____",names(pheno),value=T)
  # presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
  data_pa <- decostand(x = pheno[,species.names,with=F], "pa")
  # calculate sum per species
  data_sum <- data.frame(prevalence=apply(data_pa, 2, sum)/nrow(data_pa))
  data_sum$MGS = rownames(data_sum)
  
  a = data_sum$MGS[data_sum$prevalence <= 1/100] # 1,602 species have a prevalence greater than 1%
  pheno <- pheno[ , -a, with=F] 

# Calculate Shannon diversity ####

  mgs.names = grep("____",names(pheno),value=T) # vector with MGS names 
  pheno[,shannon:=diversity(pheno[,mgs.names, with=F],index="shannon")]
  
# Save data

  saveRDS(pheno, file=paste0(input,"pheno.MGS.Upp.rds"))
  
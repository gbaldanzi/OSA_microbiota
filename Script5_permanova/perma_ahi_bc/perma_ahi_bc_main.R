# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: April 2021

# Last update: Jan 25, 2022

# This code will investigate the beta-diversity (Bray Curtis Dissimilarity) 
# in relation to OSA severity groups in 3 different models. 
# Analysis are run using a PERMANOVA approach


# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)

output = "/home/baldanzi/Sleep_apnea/Results/"

  # Importing data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

# Importing BC matrix 
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv', header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)

source('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/permanova.fun.R')

# Transforming two-level factor variables into numeric variables 
dades = copy(pheno[valid.ahi=='yes',])
a= c("Sex","ppi","metformin","hypermed","dyslipmed")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Making sure that BC and dades have the same order of observations 
  dades <-  dades[match(rownames(BC),dades$SCAPISid),]

# Outcome - character name (length=1) with matrix distance 
  outc = "BC"

# Main Exposure - character name (length=1)
  expo = "OSAcat"


# Runing PERMANOVA in parallel ####
  source('Script5_permanova/perma_basic.model.R')

  fwrite(res, file = paste0(output,"permanova_basic.model_ahi_bc.tsv"), sep="\t")

  source('Script5_permanova/perma_full.model.R')

  fwrite(res, file = paste0(output,"permanova_full.model_ahi_bc.tsv"), sep="\t")
  
  source('Script5_permanova/perma_full.model.BMI.R')

  fwrite(res, file = paste0(output,"permanova_full.model.bmi_ahi_bc.tsv"), sep="\t")
  


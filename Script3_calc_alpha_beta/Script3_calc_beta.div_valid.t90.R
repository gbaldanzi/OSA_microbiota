# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics continued:

# This script will calculate beta diversity and on individuals with valid T90 measuremeny

rm(list=ls())

# Loading packages
library(ape)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")


# Import data
  valid.t90 <- readRDS("valid.t90_MGS.shannon_Upp.rds")

# Bray-curtis dissimilarity #### 
  
# Create dataset containing exclusively the MGS as columns/variables 
  a=grep("____",names(valid.t90),value=T)
  MGS=valid.t90[,a, with=F]
  
# Estimating the BC index - creates a matrix with the BC index between individual samples 
  # MGS as relative abundances
  print("calculating BC dissimilarity index - this will take a while")
  BC=as.matrix(vegdist(MGS,method="bray"))
  rownames(BC)=colnames(BC)=valid.t90$SCAPISid
  fwrite(BC,'/home/baldanzi/Datasets/sleep_SCAPIS/T90.BCmatrix.csv',sep=",")

# PCoA of BC
  print("running the PCoA for BC index in participants with valid T90 - this will take a while")
  pcoa.bray.t90=pcoa(BC) 
  save(pcoa.bray.t90, file = 'pc_BC_t90')

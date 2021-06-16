# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics continued:

# This script will calculate beta diversity and on individuals with valid ahi

rm(list=ls())

# Loading packages
library(ape)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")

# Import data
valid.ahi <- readRDS("validsleep_MGS.shannon_Upp.rds")

# Bray-curtis dissimilarity #### 

# Create dataset containing exclusively the MGS as columns/variables 
a=grep("____",names(valid.ahi),value=T)
MGS=valid.ahi[,a, with=F]

# Estimating the BC index - creates a matrix with the BC index between individual samples 
# MGS as relative abundances
BC=as.matrix(vegdist(MGS,method="bray"))
rownames(BC)=colnames(BC)=valid.ahi$SCAPISid
fwrite(BC,'/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv',sep=",")

#BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv',sep=",")
#BC = as.matrix(BC)
#rownames(BC) = colnames(BC) 

# Principal coordinates analysis based on BC 
pcoa.bray <- pcoa(BC)
save(pcoa.bray, file = 'pc_BC')

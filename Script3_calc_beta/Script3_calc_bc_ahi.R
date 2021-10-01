# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Latest update: 20201-09-23

# This script will calculate beta diversity and on individuals with valid ahi
  message("Script3_calc_bc_ahi.R") 

# Bray-curtis dissimilarity #### 

# Create dataset containing exclusively the MGS as columns/variables 
  mgs.names=grep("____",names(pheno),value=T)
  MGSmatrix=pheno[valid.ahi=='yes',mgs.names, with=F]

# Estimating the BC index - creates a matrix with the BC index between individual samples 
# MGS as relative abundances
BC=as.matrix(vegdist(MGSmatrix,method="bray"))
rownames(BC)=colnames(BC)=pheno[valid.ahi=='yes',SCAPISid]
fwrite(BC,'/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv',sep=",")

#BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv',sep=",")
#BC = as.matrix(BC)
#rownames(BC) = colnames(BC) 

# Principal coordinates analysis based on BC 
pcoa.bray <- pcoa(BC)
save(pcoa.bray, file = '/home/baldanzi/Datasets/sleep_SCAPIS/pc_BC')

# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Latest update: 2022-02-01

# This script will calculate beta diversity for individuals with valid AHI

# Bray-curtis dissimilarity #### 

# Create dataset containing exclusively the MGS as columns/variables 
  species.names=grep("____",names(pheno),value=T)
  MGSmatrix=pheno[valid.ahi=='yes',species.names, with=F]
  rownames(MGSmatrix) <- pheno[valid.ahi=="yes", SCAPISid]

# Estimating the BC index - creates a matrix with the BC index between individual samples 
# MGS as relative abundances
  BC <-  as.matrix(vegdist(MGSmatrix,method="bray"))
  
  rownames(BC) <- colnames(BC) <- pheno[valid.ahi=='yes',SCAPISid]
  
  fwrite(BC,file = paste0(input, 'OSA.BCmatrix.csv') ,sep=",")


  # Principal coordinates analysis based on BC 
  pcoa.bray <- pcoa(BC)
  save(pcoa.bray, file = paste0(input,'pc_BC'))

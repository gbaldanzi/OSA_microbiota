# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-01-02

# This script will calculate beta diversity and on individuals with valid T90 measurement


# Bray-curtis dissimilarity #### 
  
# Create dataset containing exclusively the MGS as columns/variables 
  mgs.names <- grep("____",names(pheno),value=T)
  MGS <- pheno[valid.t90=="yes" , mgs.names, with=F]
  
# Estimating the BC index - creates a matrix with the BC index between individual samples 
  # MGS as relative abundances
  BC <- as.matrix(vegdist(MGS, method = "bray"))
  rownames(BC) <-  colnames(BC) <- pheno[valid.t90=="yes", SCAPISid]
  
  fwrite(BC, file = paste0(input,'T90.BCmatrix.csv'), sep=",")

# PCoA of BC
  pcoa.bray.t90=pcoa(BC) 
  save(pcoa.bray.t90, file = paste0(input,'pc_BC_t90'))

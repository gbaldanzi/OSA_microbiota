# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-09-21

# Last update: 2021-09-21

# Script to run PCA on MGS. 
  library(data.table)

  # Import data 
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  
  # Matrix of MGS   
  mgs=grep("____",names(pheno),value=T)
  
  mgs.relative.abundance_valid.ahi <- pheno[valid.ahi=='yes',mgs,with=F]
  
  mgs.relative.abundance_valid.t90 <- pheno[valid.t90=='yes',mgs,with=F]
  
  # PCA
  set.seed(1)
  mgs.pca_valid.ahi <- prcomp(mgs.relative.abundance_valid.ahi)
  
  saveRDS(mgs.pca_valid.ahi, file = "/home/baldanzi/Datasets/sleep_SCAPIS/mgs.pca_valid.ahi.rds")
  
  set.seed(1)
  mgs.pca_valid.t90 <- prcomp(mgs.relative.abundance_valid.t90)
  
  saveRDS(mgs.pca_valid.t90, file = "/home/baldanzi/Datasets/sleep_SCAPIS/mgs.pca_valid.t90.rds")
  
  
  
  
  
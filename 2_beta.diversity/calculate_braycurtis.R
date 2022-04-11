# Calculate Beta-diversity (Bray-curtis dissimilarity)

# This script produces Bray-curtis dissimilary matrix and perform PCoA

  library(data.table)
  library(vegan)
  library(ape)

    
  # # Defining folders
  input <- "/home/baldanzi/Datasets/sleep_SCAPIS/"

  # Import data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
   
  # Create dataset containing exclusively the MGS as columns/variables 
  species.names <- grep("____",names(pheno),value=T)
  
  # Restricting the species matrix to participants with valid AHI and valid T90/ODI data
  MGSmatrix <- pheno[valid.ahi=="yes" | valid.t90 =="yes", species.names, with=F]
  rownames(MGSmatrix) <- pheno[valid.ahi=="yes" | valid.t90 =="yes", SCAPISid]
  
  # Create the Bray-Curtis matrix 
  
  BC <-  as.matrix(vegdist(MGSmatrix,method="bray")) # MGS as relative abundances
  
  rownames(BC) <- colnames(BC) <- rownames(MGSmatrix)
  
  
  
  # Principal coordinates analysis (PCoA)
  
  # Participants with valid AHI data 
  pcts <- pheno[valid.ahi=="yes",SCAPISid]
  pcoa.bray <- pcoa(BC[pcts, pcts])
  saveRDS(pcoa.bray, file = paste0(input,'pcoa_BC_AHI.rds'))
  
  # Participants with valid T90/ODI data 
  pcts <- pheno[valid.t90=="yes",SCAPISid]
  pcoa.bray <- pcoa(BC[pcts, pcts])
  saveRDS(pcoa.bray, file = paste0(input,'pcoa_BC_T90.ODI.rds'))
  
  # Save the BC matrix 
  BC.rownames <- rownames(BC)
    
  BC <- data.frame(BC)
  setDT(BC)
  
  BC[,rownames := BC.rownames]
  
  fwrite(BC, file = paste0(input, 'BCmatrix.tsv') )

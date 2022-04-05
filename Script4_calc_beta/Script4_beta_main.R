# Calculate Beta-diversity (Bray-curtis dissimilarity)

# Run this script to calculate Bray-curtis dissimilary matrix and produce PCoA plots

  library(data.table)
  library(vegan)
  library(ape)

    
  # Folder where dataset is located 
  input <- "/home/baldanzi/Datasets/sleep_SCAPIS/"

  # Import data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

  # Calculate Bray-curtis 
   
  # Create dataset containing exclusively the MGS as columns/variables 
  species.names=grep("____",names(pheno),value=T)
  MGSmatrix <- pheno[valid.ahi=="yes" | valid.t90 =="yes", species.names, with=F]
  rownames(MGSmatrix) <- pheno[valid.ahi=="yes" | valid.t90 =="yes", SCAPISid]
  
  # Estimating the BC index - creates a matrix with the BC index between individual samples 
  # MGS as relative abundances
  BC <-  as.matrix(vegdist(MGSmatrix,method="bray"))
  
  rownames(BC) <- colnames(BC) <- rownames(MGSmatrix)
  
  #write.csv(BC,file = paste0(input, 'BCmatrix.csv'), row.names = T)
  
  # Principal coordinates analysis 
  
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
    
  setDT(BC)
  
  BC[,rownames := BC.rownames]
  
  fwrite(BC, file = paste0(input, 'BCmatrix.tsv') )

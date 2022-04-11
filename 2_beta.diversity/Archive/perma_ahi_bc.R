# This code will investigate the beta-diversity (Bray Curtis Dissimilarity) 
# in relation to AHI severity groups 

# If 16 nodes are used to run this script, it will take around 36h. 

# Importing BC matrix 
  BC = fread(paste0(input,'BCmatrix.tsv'))
  BCrownames <- BC$rownames
  BC$rownames <- NULL
  
  BC <-  as.matrix(BC)
  colnames(BC) <- rownames(BC) <- BC_rownames

# Main Exposure - character name (length=1)
  expo = "OSAcat"

  
# Runing PERMANOVA in parallel ####
  
  # Main model 
  res <- Permanova.parallel.FUN(y = BC, exposure=expo, 
                                data = pheno[valid.ahi=="yes",], 
                                model = basic.model, nod=16)

  fwrite(res, file = paste0(results.folder,"permanova_main.model_ahi_bc.tsv"))

  
  # Extended model
  res <- Permanova.parallel.FUN(y = BC, exposure=expo, 
                                data = pheno[valid.ahi=="yes",], 
                                model = extended.model, nod=16)
  
  fwrite(res, file = paste0(results.folder,"permanova_extended.model_ahi_bc.tsv"))
  


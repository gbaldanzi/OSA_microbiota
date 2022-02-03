# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: Jan 25, 2022

# This code will investigate the beta-diversity (Bray Curtis Dissimilarity) 
# in relation to OSA severity groups in 3 different models. 
# Analysis are run using PERMANOVA 

# If 16 nodes are used to run this script, it will take around 36h. 

# Importing BC matrix 
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv', header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)

  # Data set 
  dades <- pheno[valid.ahi=="yes",]
  
# Making sure that BC and dades have the same order of observations 
  dades <-  dades[match(rownames(BC),dades$SCAPISid),]

# Main Exposure - character name (length=1)
  expo = "OSAcat"

  
# Runing PERMANOVA in parallel ####
  
  # Basic model 
  res <- Permanova.parallel.FUN(outcome = "BC", exposure=expo, 
                                data = dades, model = basic.model, nod=16)

  fwrite(res, file = paste0(output,"permanova_basic.model_ahi_bc.tsv"), sep="\t")

  
  # Full model
  res <- Permanova.parallel.FUN(outcome = "BC", exposure=expo, 
                                data = dades, model = full.model, nod=16)
  
  fwrite(res, file = paste0(output,"permanova_full.model_ahi_bc.tsv"), sep="\t")
  
  
  # Full model + BMI
  res <- Permanova.parallel.FUN(outcome = "BC", exposure=expo, 
                                data = dades, model = full.model_BMI, nod=16)
  
  fwrite(res, file = paste0(output,"permanova_full.model.bmi_ahi_bc.tsv"), sep="\t")
  


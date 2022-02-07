# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last Update: 2022-02-02

# This code will investigate the beta-diversity (Bray Curtis Dissimilarity)
# in relation to ODI in 3 different models. 
# Analysis are run using  PERMANOVA 

  # Importing BC matrix 
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/T90.BCmatrix.csv', header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)

  # Data set 
  dades <- pheno[valid.t90=="yes",]

  # Making sure that BC and dades have the same order of observations 
  dades = dades[match(rownames(BC),dades$SCAPISid),]

# Outcome - character name (length=1) with matrix distance 
outc = "BC"

# Main Exposure - character name (length=1)
expo = "odicat"

  
# Runing PERMANOVA in parallel ####

  # Basic model 
  res <- Permanova.parallel.FUN(outcome = "BC", exposure=expo, 
                              data = dades, model = basic.model, nod=16)
  
  fwrite(res, file = paste0(results.folder,"permanova_basic.model_odi_bc.tsv"), sep="\t")
  
  
  # Full model
  res <- Permanova.parallel.FUN(outcome = "BC", exposure=expo, 
                                data = dades, model = full.model, nod=16)
  
  fwrite(res, file = paste0(results.folder,"permanova_full.model_odi_bc.tsv"), sep="\t")


 
  
  

# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: April 2021
# Update: 2021-12-09 

# This code will investigate the beta-diversity (Bray Curtis Dissimilarity)
# in relation to T90% severity categories 
# Analysis are run using a PERMANOVA approach

#parallel nodes (1 for each model). Therefore, this code requires at least 3 nodes. 

# Results saved at the folder: "/home/baldanzi/Sleep_apnea/Results/"
# File model 3 - permanova_model3_t90_bc.tsv


# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)

  # Importing data

  pheno <- readRDS(file="/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")


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
  expo = "t90cat"

  
  # Runing PERMANOVA in parallel ####
  
  # Basic model 
  # res <- Permanova.parallel.FUN(outcome = "BC", exposure=expo, 
  #                              data = dades, model = basic.model, nod=16)
  
  # fwrite(res, file = paste0(results.folder,"permanova_basic.model_t90_bc.tsv"), sep="\t")
  
  
  # Full model
   res <- Permanova.parallel.FUN(outcome = "BC", exposure=expo, 
                                data = dades, model = full.model, nod=16)
  
  fwrite(res, file = paste0(results.folder,"permanova_full.model_t90_bc.tsv"), sep="\t")
  


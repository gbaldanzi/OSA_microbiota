# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: 2022-01-19
# Last Update: - 2022-01-19

# This code will run pairwise comparisons between ODI severity group in 
# relation to the beta-diversity (Bray Curtis Dissimilarity) 

# Analysis are run using a PERMANOVA approach


  # Importing data

  dades <- copy(pheno[valid.t90=="yes",])


  # Importing BC matrix 
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/T90.BCmatrix.csv', header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)

  source('Script0_functions/perma.pairwise.fun.R')


# Making sure that BC and dataset have the same order of observations 
  dades <-  dades[match(rownames(BC),dades$SCAPISid),]
  

# Runing PERMANOVA in parallel ####
  
  list.res <- pairwise.perma.fun(outcome="BC", group_var="odicat", covari=full.model, data=dades, nodes=16)
  
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results_odi.rds"))
#---------------------------------------------------------------------------#

  # Produce a final summary results table with FDR-p-values 
  
  list.res = readRDS("/home/baldanzi/Sleep_apnea/Results/pairwise.perma.results_odi.rds")
  
  list.res$table <- clean.res(list.res[1:6]) 
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results_odi.rds"))
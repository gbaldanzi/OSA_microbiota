# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Version 1: 2022-01-19
# Last Update: - 2022-01-19

# This code will investigate pairwise comparison between groups of different OSA severity based on AHI


# Loading packages 
  pacman::p_load(data.table, vegan,parallel)


  output = "/home/baldanzi/Sleep_apnea/Results/"
  

  # Importing data
  dades <- copy(pheno[valid.ahi=="yes",])


  # Importing BC matrix 
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv', header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)

  source('Script0_functions/perma.pairwise.fun.R')


# Making sure that BC and dataset have the same order of observations 
  dades <-  dades[match(rownames(BC),dades$SCAPISid),]
  
  dades[OSAcat=="no OSA", OSAcat:="noOSA"]
 

  # Runing PERMANOVA in parallel ####
  
  list.res <- pairwise.perma.fun(outcome="BC", group_var="OSAcat", covari=full.model, data=dades, nodes=16)
  

  saveRDS(list.res,file=paste0(output,"pairwise.perma.results_ahi.rds"))
#---------------------------------------------------------------------------#

  # Produce a final summary results table with FDR-p-values 
  
  list.res = readRDS(paste0(output,"pairwise.perma.results_ahi.rds"))
  
  list.res$table <- clean.res(list.res[1:6]) 
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results_ahi.rds"))
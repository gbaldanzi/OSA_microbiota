# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last Update: - 2022-02-02

# This code will investigate pairwise comparison between groups 
# of different OSA severity based on AHI

# Preparation 
source('Script5_permanova/pre_pairwise.R')

  # Importing data
  dades <- copy(pheno[valid.ahi=="yes",])


  # Importing BC matrix 
  BC = fread(paste0(input,'OSA.BCmatrix.csv'), header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)


# Making sure that BC and dataset have the same order of observations 
  dades <-  dades[match(rownames(BC),dades$SCAPISid),]
  

  # Running PERMANOVA in parallel ####
  
  dades[OSAcat=="no OSA", OSAcat:="noOSA"]
  dades[,OSAcat := factor(OSAcat, levels = c("noOSA", "Mild", "Moderate", "Severe"), 
                          labels = c("noOSA", "Mild", "Moderate", "Severe"))]
  
  list.res <- pairwise.perma.fun(outcome="BC", group_var="OSAcat", covari=full.model, data=dades, nodes=16)
  

  saveRDS(list.res,file=paste0(output,"pairwise.perma.results_ahi.rds"))
#---------------------------------------------------------------------------#

  # Produce a final summary results table with FDR-p-values 
  
  list.res = readRDS(paste0(output,"pairwise.perma.results_ahi.rds"))
  
  list.res$table <- clean.res(list.res[1:6]) 
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results_ahi.rds"))
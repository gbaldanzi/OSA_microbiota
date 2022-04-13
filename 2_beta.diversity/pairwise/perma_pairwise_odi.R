# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# Version 1: 2022-01-19
# Last Update: - 2022-02-02

# This code will run pairwise comparisons between ODI severity group in 
# relation to the beta-diversity (Bray Curtis Dissimilarity) 

# Preparation 
source('2_beta.diversity/pairwise/pre_pairwise.R')

  # Importing data
  dades <- copy(pheno[valid.t90=="yes",])

# Making sure that BC and dataset have the same observations and same order
  BC <- BC[dades$SCAPISid,dades$SCAPISid]
  
  dades <-  dades[match(rownames(BC),dades$SCAPISid),]
  

# Runing PERMANOVA in parallel ####
  
  list.res <- pairwise.perma.fun(outcome="BC", group_var="odicat", covari=main.model, data=dades, nodes=16)
  
  
  saveRDS(list.res,file=paste0(results.folder,"pairwise.perma.results_odi.rds"))
#---------------------------------------------------------------------------#

  # Produce a final summary results table with FDR-p-values 
  
  list.res = readRDS(paste0(results.folder,"pairwise.perma.results_odi.rds"))
  
  list.res$table <- clean.res(list.res[1:6]) 
  
  saveRDS(list.res,file=paste0(results.folder,"pairwise.perma.results_odi.rds"))
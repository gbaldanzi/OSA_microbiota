# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# Last Update: - 2022-02-02

# This code will run pairwise comparisons between T90 severity group in 
# relation to the beta-diversity (Bray Curtis Dissimilarity) 

# Preparation 
source('2_beta.diversity/pairwise/pre_pairwise.R')

  # Importing data
    dades <- copy(pheno[valid.t90=="yes",])


  # Making sure that BC and dataset have the same order of observations 
    BC <- BC[dades$SCAPISid,dades$SCAPISid]
    
    dades <-  dades[match(rownames(BC),dades$SCAPISid),]

  # Running PERMANOVA in parallel ####
  
  dades[, t90cat:=factor(t90cat, levels(t90cat), 
                         labels = c("t0", "t1", "t2", "t3"))]
  
  
  list.res <- pairwise.perma.fun(outcome="BC", group_var="t90cat", covari=main.model, data=dades, nodes=16)
  
  
  saveRDS(list.res,file=paste0(results.folder,"pairwise.perma.results_t90.rds"))
#---------------------------------------------------------------------------#

  # Produce a final summary results table with FDR-p-values 
  
  list.res = readRDS(paste0(results.folder,"pairwise.perma.results_t90.rds"))
  
  list.res$table <- clean.res(list.res[1:6]) 
  
  saveRDS(list.res,file=paste0(results.folder,"pairwise.perma.results_t90.rds"))
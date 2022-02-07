# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# Last Update: - 2022-02-02

# This code will run pairwise comparisons between T90 severity group in 
# relation to the beta-diversity (Bray Curtis Dissimilarity) 

# Preparation 
source('Script5_permanova/pre_pairwise.R')

  # Importing data
    dades <- copy(pheno[valid.t90=="yes",])

  # Importing BC matrix 
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/T90.BCmatrix.csv', header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)


# Making sure that BC and dataset have the same order of observations 
  dades <-  dades[match(rownames(BC),dades$SCAPISid),]
  

# Runing PERMANOVA in parallel ####
  
  dades[, t90cat:=factor(t90cat, levels(t90cat), 
                         labels = c("t0", "t1", "t2", "t3"))]
  
  
  list.res <- pairwise.perma.fun(outcome="BC", group_var="t90cat", covari=full.model, data=dades, nodes=16)
  
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results_t90.rds"))
#---------------------------------------------------------------------------#

  # Produce a final summary results table with FDR-p-values 
  
  list.res = readRDS("/home/baldanzi/Sleep_apnea/Results/pairwise.perma.results_t90.rds")
  
  list.res$table <- clean.res(list.res[1:6]) 
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results_t90.rds"))
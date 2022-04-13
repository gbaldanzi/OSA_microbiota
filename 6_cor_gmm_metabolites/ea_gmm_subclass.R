# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script will run the enrichment analysis for metabolites groups (Sub-pathways) 
# in the metabolites associations with T90-associated GMM. 

  library(data.table)
  library(fgsea)
  library(tidyverse)
  library(BiocParallel)
  library(rio)

  results.folder = "/home/baldanzi/Sleep_apnea/Results/"
  input <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/data_processed/'
  
  set.seed(10)
  cores <- 16
  
  # gmm-metabolties correlation results 
  correlations <- import(paste0(results.folder,"cor_gmm_metabolites.tsv"))
  correlations <- correlations[which(!is.na(correlations$p.value)), ]
  
  # subclass list 
  subclass <-  readRDS(paste0(input,'subclass.list.rds'))
  
  
  # enrichment function for metabolite subclasses
  
  subclass.fun <- function(data) {
    
    res <- bplapply(unique(data$modules), function(modules) {
      
      # rank based on p-value
      
      data <- data[which(data$modules == modules), ]
      stats <- rank(-data$p.value, na = "keep")
      names(stats) <- data$metabolites
      
      # run gsea
      
      res <- as.data.frame(fgsea(subclass, stats, scoreType = "pos", eps = 0))
      data.frame(modules = modules, subclass = res$pathway, estimate = res$NES, p.value = res$pval, 
                 q.value = res$padj,
                 size = res$size, leading = sapply(res$leadingEdge, function(x) paste(x, collapse = ";")))
      
    }, BPPARAM = MulticoreParam(cores))
    
    do.call(rbind, res)
    
  }
  
  # run gsea for postive correlations 
  subclass.pos <- subclass.fun(correlations[which(correlations$rho >= 0), ])
  subclass.pos$direction <- "positive"
  
  # run gsea for negative correlations
  subclass.neg <- subclass.fun(correlations[which(correlations$rho < 0), ])
  subclass.neg$direction <- "negative"
  
  # combine results 
  
  subclass.res <- rbind(subclass.pos, subclass.neg)
  
  # save results
  fwrite(subclass.res[,c("modules", "subclass", "estimate", "p.value", "q.value","leading", "size","direction")], 
         file = paste0(results.folder, "ea_gmm_subclass.tsv"))
  
  
  
  
  
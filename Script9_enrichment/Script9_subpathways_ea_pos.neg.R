# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-09-16

# Last update: 2021-09-16

# Enrichment analysis of sub-pathways among associations between 
# signature species and the metabolites

# MGS-metabolites associations stratified by the direction of cor. coef. 

# Correlation between relevant MGS and metabolites 
pacman::p_load(data.table,fgsea,tidyverse)

rm(list = ls())

input = "/home/baldanzi/Sleep_apnea/Results/"
output = "/home/baldanzi/Sleep_apnea/Results/"

  # Import results 
  mgs_met <- fread("/home/baldanzi/Datasets/Mgs_metab_correlations/table1_shortmgs.tsv") #Correlation metabolites and MGS

  res.m2 <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  mgs.fdr.m2 <- readRDS('/home/baldanzi/Sleep_apnea/Results/mgs.m2.rds')
  mgs.fdr.m2 <- unique(do.call('c',mgs.fdr.m2))   #Maintax_mgs
  
  # MGS correlated to either AHI, ODI or T90 but not to BMI (relevant mgs)
  mgs.rel <- res.m2 %>% filter(exposure =="ahi" | exposure == "t90" | exposure =="odi") %>%
              filter(MGS %in% mgs.fdr.m2) %>% select(mgs)
  mgs.rel <- unique(mgs.rel$mgs) # 42 mgs (HG3A)
  
  # Restricting MGS-metabolites correlation results 
  
  sum(mgs.rel %in% mgs_met[,mgs]) # 37 MGS are included in the Atlas
  
  mgs.rel_met <- mgs_met[mgs %in% mgs.rel,]
  
  mgs.rel_met_pos <- mgs.rel_met[correlation.sp>0,] 
  
  mgs.rel_met_neg <- mgs.rel_met[correlation.sp<0,] 
  
  mgs.rel[!mgs.rel %in% mgs_met[,mgs]] # 5 mgs used in the sleep apnea study were not used
  # in the metabolomics study 
  
  m <- mgs.rel[mgs.rel %in% mgs_met[,mgs]]
  
  metabolites.pvalues_pos <- lapply(m, function(x) {
    
    p <- mgs.rel_met_pos[mgs == x,p.value.sp]
    names(p) <- mgs.rel_met_pos[mgs == x,metabolite]
    return(p)
    
    })
  
  metabolites.pvalues_neg <- lapply(m, function(x) {
    
    p <- mgs.rel_met_neg[mgs == x,p.value.sp]
    names(p) <- mgs.rel_met_neg[mgs == x,metabolite]
    return(p)
    
  })
  
    names(metabolites.pvalues_pos) <- m
    names(metabolites.pvalues_neg) <- m
    
    # Enrichment analysis 
    
    subpathways <-  readRDS('/home/baldanzi/Datasets/Mgs_metab_correlations/subpathwayslist.rds')
    
    # Positive correlations
    set.seed(123)
    res_pos <- lapply(m,function(x){
                  temp <- fgsea(pathways = subpathways, stats = rank(-metabolites.pvalues_pos[[x]]), 
                        eps=0,scoreType="pos", nPermSimple = 10000, nproc=4)
                  temp[,mgs:=x]
    })
    
    names(res_pos) <- m
    
    
    res.table_pos <- do.call(rbind,res_pos)
    
    # Negative correlations
    set.seed(123)
    res_neg <- lapply(m,function(x){
      temp <- fgsea(pathways = subpathways, stats = rank(-metabolites.pvalues_neg[[x]]), 
                    eps=0,scoreType="pos", nPermSimple = 10000, nproc=4)
      temp[,mgs:=x]
    })
    
    names(res_neg) <- m
    
    
    res.table_neg <- do.call(rbind,res_neg)
    
    # Adjusting names to fullnames 
    
    fullnames <- unique(res.m2[mgs %in% m,.(mgs,MGS)])
    
    res.table_pos <- merge(res.table_pos,fullnames , by="mgs", all.x=T, all.y=F)
    res.table_neg <- merge(res.table_neg,fullnames , by="mgs", all.x=T, all.y=F)
    
    # Save results 
    
    fwrite(res.table_pos, file= paste0(output,"ea_subpathways_pos.tsv"), sep='\t')
    fwrite(res.table_neg, file= paste0(output,"ea_subpathways_neg.tsv"), sep='\t')
  
  
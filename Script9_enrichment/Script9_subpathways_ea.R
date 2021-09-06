# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-13

# Last update: 2021-09-03

# Enrichment analysis of subpathways among the metabolites correlated to the identified MGSs

# Correlation between relevant MGS and metabolites 
pacman::p_load(data.table,fgsea,tidyverse)

rm(list = ls())

input = "/home/baldanzi/Sleep_apnea/Results/"
output = "/home/baldanzi/Sleep_apnea/Results/"

  # Import results 
  mgs_met <- fread("/home/baldanzi/Datasets/Mgs_metab_correlations/table1_shortmgs.tsv") #Correlation metabolites and MGS

  res.m2 <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  
  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]

  # Select the relevant MGSs 
  mgs.bmi <- res.m2[exposure =="BMI" & q.value<.05,mgs] # MGS correlated to BMI
  
  # MGS correlated to either AHI or T90 but not to BMI (relevant mgs)
  mgs.rel <- res.m2 %>% filter(exposure =="ahi" | exposure == "t90") %>% 
    filter(q.value<.05) %>% filter(!mgs %in% mgs.bmi) %>% select(mgs)
  mgs.rel <- unique(mgs.rel$mgs)
  
  # Restricting MGS-metabolites correlation results 
  
  mgs.rel_met <- mgs_met[mgs %in% mgs.rel,]
  
  mgs.rel[!mgs.rel %in% mgs_met[,mgs]] # 6 mgs used in the sleep apnea study were not used
  # in the metabolomics study 
  
  m <- mgs.rel[mgs.rel %in% mgs_met[,mgs]]
  
  metabolites.pvalues <- lapply(m, function(x) {
    
    p <- mgs.rel_met[mgs == x,p.value.sp]
    names(p) <- mgs.rel_met[mgs == x,metabolite]
    return(p)
    
    })
  
    names(metabolites.pvalues) <- m
    
    # Enrichment analysis 
    
    subpathways <-  readRDS('/home/baldanzi/Datasets/Mgs_metab_correlations/subpathwayslist.rds')
    
    
    set.seed(123)
    res <- lapply(m,function(x){
                  temp <- fgsea(pathways = subpathways, stats = rank(-metabolites.pvalues[[x]]), 
                        eps=0,scoreType="pos", nPermSimple = 10000, nproc=4)
                  temp[,mgs:=x]
    })
    
    names(res) <- m
    
    
    res.table <- do.call(rbind,res)
    
    fullnames <- unique(res.m2[mgs %in% m,.(mgs,MGS)])
    
    res.table <- merge(res.table,fullnames , by="mgs", all.x=T, all.y=F)
    
    fwrite(res.table, file= paste0(output,"ea_subpathways.tsv"), sep='\t')
    
  
  
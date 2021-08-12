# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-05

# This script runs the overrepresentation analysis for metabolites sub_pathways among the metabolites
# correlated to MGS identified in model 1 

  library(fgsea)
  library(data.table)

  subpathways <-  readRDS('/home/baldanzi/Datasets/Mgs_metab_correlations/subpathwayslist.rds')
  list_mm_p.value.sp <-  readRDS('/home/baldanzi/Sleep_apnea/Results/list_mm_p.value.sp.rds')
  
  set.seed(123)
  res.ahi <- fgsea(pathways = subpathways, stats = rank(-list_mm_p.value.sp[[1]]), eps=0,scoreType="pos", nPermSimple = 100000, nproc=3)
  res.ahi=cbind("AHI",res.ahi)
  saveRDS(res.ahi,'/home/baldanzi/Sleep_apnea/Results/temp_subpath_ea_ahi.rds')
  
  res.bmi <- fgsea(pathways = subpathways, stats = rank(-list_mm_p.value.sp[[2]]), eps=0,scoreType="pos", nPermSimple = 100000)
  res.bmi=cbind("BMI",res.bmi)
  saveRDS(res.bmi,'/home/baldanzi/Sleep_apnea/Results/temp_subpath_ea_bmi.rds')
  
  res.t90 <- fgsea(pathways = subpathways, stats = rank(-list_mm_p.value.sp[[3]]), eps=0,scoreType="pos", nPermSimple = 100000)
  res.t90=cbind("T90",res.t90)
  saveRDS(res.t90,'/home/baldanzi/Sleep_apnea/Results/temp_subpath_ea_t90.rds')
  
  res.subpath.ea = rbind(res.ahi,res.bmi,res.t90)
  names(res.subpath.ea)[1] <- "exposure"
  setnames(res.subpath.ea,"padj","q.value")
  fwrite(res.subpath.ea, '/home/baldanzi/Sleep_apnea/Results/ea_subpath.tsv',sep='\t')
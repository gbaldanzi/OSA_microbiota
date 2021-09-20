# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-02

# Last update: 2021-09-03

# This script imports the table of correlations from the paper of Koen et al. and 
# filters to the MGS identified with model 1 of this project. Finally, it saves the sp 
# p.value corresponding to MGS in a list file. ("/home/baldanzi/Sleep_apnea/Results/list_mm_p.value.sp.rds")

library(data.table)
 #mgs_met <- fread("/home/baldanzi/Datasets/Mgs_metab_correlations/table1.tsv")
 #mgs_met[,mgs:=substring(metagenomic.species,1,9)]
 #setcolorder(mgs_met,"mgs")
 #fwrite(mgs_met,"/home/baldanzi/Datasets/Mgs_metab_correlations/table1_shortmgs.tsv")
mgs_met <- fread("/home/baldanzi/Datasets/Mgs_metab_correlations/table1_shortmgs.tsv")

# MGSs identified in model 1

# Importing results 
input = "/home/baldanzi/Sleep_apnea/Results/"
output = "/home/baldanzi/Sleep_apnea/Results/"

  res <- fread(paste0(input,"cor_all.var_mgs.tsv"))
  res[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  res.list = list(AHI = res[exposure=="ahi",],
                  T90 = res[exposure=="t90",],
                  BMI = res[exposure=="BMI",])

  # filter MGS significant at the FDR p-value<0.05
  mgs.fdr = lapply(res.list,function(x){x[q.value<0.05,mgs]})


  # List mgs_meta correlations p.value and metabolites names by phenotype 
  list_mm_p.value.sp <- list(AHI=NULL, T90=NULL, BMI=NULL)
  for(i in names(mgs.fdr)){
  message(i)
  p <- mgs_met[mgs %in% mgs.fdr[[i]],p.value.sp]
  list_mm_p.value.sp[[i]] <- as.numeric(p)
  names(list_mm_p.value.sp[[i]]) <- mgs_met[mgs %in% mgs.fdr[[i]],metabolite]
  }

  saveRDS(list_mm_p.value.sp,paste0(output,"list_mm_p.value.sp.rds"))


  # Pathway list 
  annotation <- fread('/home/baldanzi/Datasets/Metabolon/clean/original/scapis_annotations_clean.tsv')
  annotation[,SUB_PATHWAY:=gsub("-","_",SUB_PATHWAY)]
  annotation[,SUB_PATHWAY:=gsub("/","_",SUB_PATHWAY)]
  annotation[,SUB_PATHWAY:=gsub('(',"_",SUB_PATHWAY,fixed=T)]
  annotation[,SUB_PATHWAY:=gsub(")","_",SUB_PATHWAY)]
  annotation[,SUB_PATHWAY:=gsub(",","",SUB_PATHWAY)]
  annotation[,SUB_PATHWAY:=gsub(";","",SUB_PATHWAY)]
  annotation[,SUB_PATHWAY:=gsub(" ","_",SUB_PATHWAY)]
  lev <- annotation[,SUB_PATHWAY]
  lev <- lev[!is.na(lev)]
  lev <- lev[lev!=""]
  pathways=NULL
  pathways=list()
  for(i in lev)
  {
    eval(parse(text=paste0("pathways$", i , "=annotation[SUB_PATHWAY==i,CHEMICAL_NAME]")))
  }
  
  saveRDS(pathways, '/home/baldanzi/Datasets/Mgs_metab_correlations/subpathwayslist.rds')
  
  # Metabolites annotation  
  metabolites <- annotation[,CHEMICAL_NAME]
  metlist = as.list(metabolites)
  names(metlist) = metabolites
  
  saveRDS(metlist, '/home/baldanzi/Datasets/Mgs_metab_correlations/metabolitelist.rds')

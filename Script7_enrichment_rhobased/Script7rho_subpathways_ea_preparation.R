# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-05

# This script imports the table of correlations from the paper of Koen et al. and 
# filters to the MGS identified with model 1 of this project. Finally, it saves the sp 
# p.value corresponding to MGS in a list file. ("/home/baldanzi/Sleep_apnea/Results/list_mm_p.value.sp.rds")

library(data.table)

input = "/home/baldanzi/Sleep_apnea/Results/"
output = "/home/baldanzi/Sleep_apnea/Results/"

 #mgs_met <- fread("/home/baldanzi/Datasets/Mgs_metab_correlations/table1.tsv")
 #mgs_met[,mgs:=substring(metagenomic.species,1,9)]
 #setcolorder(mgs_met,"mgs")
 #fwrite(mgs_met,"/home/baldanzi/Datasets/Mgs_metab_correlations/table1_shortmgs.tsv")
mgs_met <- fread("/home/baldanzi/Datasets/Mgs_metab_correlations/table1_shortmgs.tsv")

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

  # MGSs identified in model 1

  res.ahi <- fread(paste0(input,"cor_ahi_mgs.tsv"))
  res.bmi <- fread(paste0(input,"cor_BMI_mgs.tsv"))
  res.t90 <- fread(paste0(input,"cor_t90_mgs.tsv"))

  res.list = list(res.ahi,res.bmi,res.t90)

  # filter MGS significant at the FDR p-value<0.05 for model 1 ####
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})

  names(mgs.fdr) <- c("AHI","BMI","T90")

  mgs.fdr = lapply(mgs.fdr,function(x){
  temp=strsplit(x, "____")
  temp=matrix(unlist(temp), ncol=2, byrow=TRUE)
  return(temp[,2])
  })
  
  # List mgs_meta correlations p.value and metabolites names by phenotype 
  list_mm_correlation.sp <- list(ahi=NULL, bmi=NULL, t90=NULL)
  for(i in 1:3){
    print(i)
    p <- mgs_met[mgs %in% mgs.fdr[[i]],correlation.sp]
    list_mm_correlation.sp[[i]] <- as.numeric(p)
    names(list_mm_correlation.sp[[i]]) <- mgs_met[mgs %in% mgs.fdr[[i]],metabolite]
  }
  
  saveRDS(list_mm_correlation.sp,paste0(output,"list_mm_correlation.sp.rds"))

  # MGSs identified in model 2 exclusively to each phenotype ####
  
  # Importing results 
  
  res.ahi <- fread(paste0(input,"cor2_ahi_mgs.tsv"))
  res.bmi <- fread(paste0(input,"cor2_BMI_mgs.tsv"))
  res.t90 <- fread(paste0(input,"cor2_t90_mgs.tsv"))
  
  res.list2 = list(res.ahi,res.bmi,res.t90)
  
  # filter MGS significant at the FDR p-value<0.05 for model 1
  mgs.fdr = lapply(res.list2,function(x){x[x$q.value<0.05,MGS]})
  
  names(mgs.fdr) <- c("AHI","BMI","T90")
  
  a <-  unlist(mgs.fdr)
  df.mgs.frd <- data.frame(exp=substring(names(a),0,3), mgs=a)
  mgs.fdr <- lapply(names(mgs.fdr), function(x){
    c <- df.mgs.fdr[df.mgs.fdr$exp!=x,"mgs"]
    return(df.mgs.fdr[df.mgs.fdr$exp==x & !df.mgs.fdr$mgs %in% c,"mgs"])
    })
  
  names(mgs.fdr) <- c("AHI","BMI","T90")
  
  mgs.fdr = lapply(mgs.fdr,function(x){
    temp=strsplit(x, "____")
    temp=matrix(unlist(temp), ncol=2, byrow=TRUE)
    return(temp[,2])
  })
  
  # List mgs_meta correlations p.value and metabolites names by phenotype 
  list2_mm_correlation.sp <- list(ahi=NULL, bmi=NULL, t90=NULL)
  for(i in 1:3){
    print(i)
    p <- mgs_met[mgs %in% mgs.fdr[[i]],correlation.sp]
    list_mm_correlation.sp[[i]] <- as.numeric(p)
    names(list_mm_correlation.sp[[i]]) <- mgs_met[mgs %in% mgs.fdr[[i]],metabolite]
  }
  
  saveRDS(list_mm_correlation.sp,paste0(output,"list2_mm_correlation.sp.rds"))


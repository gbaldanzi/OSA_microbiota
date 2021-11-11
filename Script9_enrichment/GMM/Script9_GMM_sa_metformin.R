# Project Sleep Apnea and Gut microbiota 

# Gabriel Baldanzi 

# 2021-11-09

# Last update: 2021-11-09

# The aim of this script it to conduct a sensitivity analysis for the GMM modules. 
# For the sensitivity analysis, we will also adjust for metformin given the connection between
# metformin and serum lactate. 
# In this script, we will run the correlation between AHI and MGS abudance and the correlation 
# between T90 and MGS abundance. 
# In a second step, we will investigate which GMM models are enriched in this results 


  # Loading packages 
  pacman::p_load(data.table, ppcor, fastDummies, vegan)

  # Input and output folders 
  input1 = '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/'
  output.results = "/home/baldanzi/Sleep_apnea/Results/"

  # Importing data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")


  # Removing rare MGS   
  # Calculating MGS prevalence 
  noms=grep("____",names(pheno),value=T)
  # presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
  data_pa <- decostand(x = pheno[,noms,with=F], "pa")
  # calculate sum per species
  data_sum <- data.frame(prevalence=apply(data_pa, 2, sum),
                         percentage = apply(data_pa,2,sum)/nrow(pheno))
  data_sum$MGS = rownames(data_sum)
  a = data_sum$MGS[data_sum$prevalence<5] #19 MGS are only present in less than 5 individuals
  #a = data_sum$MGS[data_sum$percentage<.01] #383 MGS are present in less than 1% of individuals
  
  pheno <- pheno[ , -a, with=F] 
  
  
  # Sensitivity analysis - remove metformin users 
  pheno <- pheno[metformin == "no",]
  
  # Correlations 
  source(paste0(input1,"Spearman.correlation.function.R")) # Correlation function 
  
  source('cor_model1/Script6_cor1_AHI_MGS.R')
  source('cor_model1/Script6_cor1_T90_MGS.R')
  
  
  res <- rbind(res.ahi, res.t90)
  
  res$model= "model1"
  
  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                 "N", "method", "covariates","q.value","model")

  
  # Enrichment analysis ####
  
  pacman::p_load(tidyr, fgsea,rio)
  
  source("MGS.Enrich.function.R")
  
  # input and output folders 
  input1 = "/home/baldanzi/Sleep_apnea/Results/"
  input2 = "/home/baldanzi/Datasets/MGS/original/"
  output = "/home/baldanzi/Sleep_apnea/Results/"

  setDT(res)
  
  cutlast <- function(char,n){
    l <- nchar(char)
    a <- l-n+1
    return(substr(char,a,l))
  }
  
  res[,mgs:=cutlast(MGS,9)]
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",])
  
  # List of modules ####
  load(paste0(input2,'MGS_HG3A.GMMs2MGS.RData')) # object = MGS_HG3A.GMMs2MGS
  list.modules <-  MGS_HG3A.GMMs2MGS
  
  # Enrichment analysis   
  # Positive correlations 
  message("Positive correlations")
  res.pos = list()[1:2]
  
  for(i in 1:2){
    print(i)
    
    res.pos[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                                         p.value.name="p.value",
                                         cor.var.name = "cor.coefficient",
                                         MGS.var.name = "mgs",
                                         enrich.var.list = list.modules,
                                         direction = "positive", 
                                         maxSize = Inf)
    names(res.pos)[i] <- names(res.list)[i]
  }
  
  
  # Merging with module annotation 
  
  res.pos <- do.call(rbind,res.pos)
  
  fwrite(res.pos, file = paste0(output,"ea_GMM_pos_sa_metformin.tsv"))
  # res.pos <- fread("/home/baldanzi/Sleep_apnea/Results/ea_GMM_pos_sa_metformin.tsv")
  
  # Clearning the results 
  
  # Import GMM names 
  gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  
  res.pos <- merge(res.pos, gmm.names, by.x="pathway", by.y="Module")

  #Cleaning
  res.pos[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  res.pos[,NES:=round(NES,3)]

  res.pos <-  res.pos[,.(exposure,pathway,Name,NES,pval,q.value)]

  #Saving 
  saveRDS(res.pos,file=paste0(output,"ea_GMM_table_pos_sa_metformin.rds"))

  
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-16

# Inferential Statistics 

# This code will investigate the association between T90 and MGS 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/cor_t90_mgs.tsv"

rm(list = ls())
# Loading packages 
pacman::p_load(data.table, ppcor, fastDummies, vegan)

# Output folders 
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
valid.t90 <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/valid.t90_MGS.shannon_Upp.rds")

#Calculating  MGS prevalence ####
noms=grep("____",names(valid.t90),value=T)
# presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
data_pa <- decostand(x = valid.t90[,noms,with=F], "pa")
# calculate sum per species
data_sum <- data.frame(prevalence=apply(data_pa, 2, sum),
                       percentage = apply(data_pa,2,sum)/nrow(valid.t90))
data_sum$MGS = rownames(data_sum)
#a = data_sum$MGS[data_sum$prevalence<5] #19 MGS are only present in less than 5 individuals
a = data_sum$MGS[data_sum$percentage<.01] # 387 MGS are present in less than 1% of individuals

# Removing MGS that are rare
valid.t90 <- valid.t90[ , -a, with=F] 

# Transforming two-level factor variables into numeric variables 
dades = copy(valid.t90)
a= c("Sex", "metformin","hypermed","dyslipmed","ppi")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming characters to factor variables 
dades[,plate:=as.factor(dades$plate)]

#Spearman correlation function ####
source("Spearman.correlation.function.R")

# Correlation between t90 and MGS - Step1 ####

#Preparing exposure, outcomes, and covariates
exposure="t90"
outcomes=grep("___",names(valid.t90),value=T)

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate + received 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")

# Running correlation 
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = model1,data = dades)

  res = res[order(res$q.value),]

  res$model= "model1"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")

  #fwrite(res, file = paste0(output,"cor_t90_mgs.tsv"), sep="\t")
  fwrite(res, file = paste0(output,"cor_t90_mgs_filter001.tsv"), sep="\t")
  
#--------------------------------------------------------------------------#
# Merging results with taxonomy information #### 

#taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
#setnames(taxonomy,"maintax_mgs","MGS")

#dades <- fread(paste0(output,"cor_t90_mgs.tsv"))
#dades <- merge(dades, taxonomy, by="MGS", all.x=T)
#fwrite(dades, file=paste0(output,"cor_t90_mgs.tsv"))
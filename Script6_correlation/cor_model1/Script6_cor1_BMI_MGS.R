# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-16

# Inferential Statistics 

# This code will investigate the association between BMI and MGS 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/cor_BMI_mgs.tsv"

rm(list = ls())
# Loading packages 
pacman::p_load(data.table, ppcor, fastDummies, vegan, ggplot2,parallel)

# Output folders 
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

#Calculate Shannon diversity ####
a = grep("____",names(pheno),value=T) # vector with MGS names 
pheno[,shannon:=diversity(pheno[,a, with=F],index="shannon")]

#Calculating  MGS prevalence ####
noms=grep("____",names(pheno),value=T)
# presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
data_pa <- decostand(x = pheno[,noms,with=F], "pa")
# calculate sum per species
data_sum <- data.frame(prevalence=apply(data_pa, 2, sum))
data_sum$MGS = rownames(data_sum)
a = data_sum$MGS[data_sum$prevalence<5] #19 MGS are only present in less than 5 individuals

# Removing MGS that are rare
pheno <- pheno[ , -a, with=F] 

# Transforming two-level factor variables into numeric variables 
dades = copy(pheno)
a= c("Sex", "metformin","hypermed","dyslipmed","ppi")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming characters to factor variables 
dades[,plate:=as.factor(dades$plate)]

#Spearman correlation function ####
source("Spearman.correlation.function.R")

# Correlation between BMI and MGS - Step1 ####

#Preparing exposure, outcomes, and covariates
exposure="BMI"
outcomes=grep("___",names(pheno),value=T)

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate + received 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")

# Running correlation 
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = model1,data = dades)

  res = res[order(res$q.value),]

  res$model= "model1"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")

fwrite(res, file = paste0(output,"cor_BMI_mgs.tsv"), sep="\t")

#--------------------------------------------------------------------------#
# Merging results with taxonomy information #### 

taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
setnames(taxonomy,"maintax_mgs","MGS")

dades <- fread(paste0(output,"cor_BMI_mgs.tsv"))
dades <- merge(dades, taxonomy, by="MGS", all.x=T)
fwrite(dades, file=paste0(output,"cor_BMI_mgs.tsv"))
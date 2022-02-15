# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script to prepare data to be used in STATA

# STATA will be used to perform analysis with the inputed data set 

# Because AHI have a lower sample size than ODI and T90, we are inputing 
# missing AHI values 

library(data.table)
library(vegan)
library(fastDummies)

# Defining output folders 
input <-  "/home/baldanzi/Datasets/sleep_SCAPIS/"

# Importing data
pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

# Removing rare species (species with prevalence < or = 1%)   
# Calculating MGS prevalence 
species.names=grep("____",names(pheno),value=T)
# presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
data_pa <- decostand(x = pheno[,species.names,with=F], "pa")
# calculate sum per species
data_sum <- data.frame(prevalence=apply(data_pa, 2, sum)/nrow(data_pa))
data_sum$MGS = rownames(data_sum)

a = data_sum$MGS[data_sum$prevalence <= 1/100] # 1,602 species have a prevalence greater than 1%
pheno <- pheno[ , -a, with=F] 

  # Create dummy variables for factor variables 

  #Covariates 
  basic.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  full.model <-  c(basic.model,"BMI", "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", 
                 "visit.month", "metformin","hypermed","dyslipmed","ppi")

  numeric.covari = names(pheno[,full.model,with=F])[sapply(pheno[,full.model,with=F],is.numeric)]
  factor.covari =  names(pheno[,full.model,with=F])[sapply(pheno[,full.model,with=F],is.factor)]
  factor.covari = c(factor.covari, names(pheno[,full.model,with=F])[sapply(pheno[,full.model,with=F],is.character)])
  
  factor.covari = names(pheno[,factor.covari,with=F])[sapply(pheno[,factor.covari,with=F],function(x) ifelse(length(unique(x[!is.na(x)]))>2,T,F))]
  
  temp.dataset = dummy_cols(pheno, select_columns = factor.covari,
                            remove_most_frequent_dummy = T, 
                            remove_selected_columns = F)
  
  setnames(temp.dataset, "leisurePA_mostly sedentary", "leisurePA_sed")
  
  dummy.names = names(temp.dataset)[!names(temp.dataset) %in% names(pheno)]
  
  a <- c("SCAPISid", dummy.names)
  temp.dataset <- temp.dataset[,a,with=F]
  
  pheno <- merge(pheno, temp.dataset, by="SCAPISid")
  

# Makes species names shorter 

cutlast <- function(char,n){
  l <- nchar(char)
  a <- l-n+1
  return(substr(char,a,l))
}


mgs.names.index <- grep("____",names(pheno))
names(pheno)[mgs.names.index] <- cutlast(names(pheno)[mgs.names.index],9)

# Exporting in STATA friendly format 



variables.export <- c("SCAPISid", "ahi", "odi", "t90", "valid.ahi",
                      full.model, dummy.names,
                      grep("HG3A",names(pheno),value=T),
                      "lowestsat")

pheno <- pheno[valid.t90=="yes",]


require(foreign)
write.dta(pheno[,variables.export,with=F], paste0(input,"pheno.dta"))
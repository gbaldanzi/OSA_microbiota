# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script to prepare data to be used in STATA

# STATA will be used to perform analysis with the inputed data set 

# Because AHI have a lower sample size than ODI and T90, we are inputing 
# missing AHI values 

rm(list=ls())

library(data.table)
library(vegan)
library(fastDummies)

# Defining output folders 
input <-  "/home/baldanzi/Datasets/sleep_SCAPIS/"

# Importing data
pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

  # Create dummy variables for factor variables 

  #Covariates 
  basic.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon","BMI")
  
  
  numeric.covari = names(pheno[,basic.model,with=F])[sapply(pheno[,basic.model,with=F],is.numeric)]
  factor.covari =  names(pheno[,basic.model,with=F])[sapply(pheno[,basic.model,with=F],is.factor)]
  factor.covari = c(factor.covari, names(pheno[,basic.model,with=F])[sapply(pheno[,basic.model,with=F],is.character)])
  
  factor.covari = names(pheno[,factor.covari,with=F])[sapply(pheno[,factor.covari,with=F],function(x) ifelse(length(unique(x[!is.na(x)]))>2,T,F))]
  
  temp.dataset = dummy_cols(pheno, select_columns = factor.covari,
                            remove_most_frequent_dummy = T, 
                            remove_selected_columns = F)
  
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
                      "visit.month",
                      basic.model, dummy.names,
                      grep("HG3A",names(pheno),value=T))

pheno <- pheno[valid.t90=="yes",]

pheno[,leisurePA:=factor(leisurePA, labels=c("PA1","PA2","PA3","PA4"))]

pheno[,educat:=factor(educat, labels=c("cat0", "cat1", "cat2", "cat3"))]


require(foreign)
write.dta(pheno[,variables.export,with=F], paste0(input,"pheno.dta"))
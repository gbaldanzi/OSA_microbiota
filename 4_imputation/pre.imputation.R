# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script to prepare data to be used in STATA

# STATA will be used to perform analysis with the imputed data set 

# Because AHI have a lower sample size than ODI and T90, we are imputing 
# missing AHI values 

  rm(list=ls())

  library(data.table)
  library(vegan)
  library(fastDummies)

  # Folders 
  work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'

  # Import data
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))

  # Create dummy variables for factor variables 

  #Covariates 
  main.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon","BMI")
  
  
  numeric.covari = names(pheno[,main.model,with=F])[sapply(pheno[,main.model,with=F],is.numeric)]
  factor.covari =  names(pheno[,main.model,with=F])[sapply(pheno[,main.model,with=F],is.factor)]
  factor.covari = c(factor.covari, names(pheno[,main.model,with=F])[sapply(pheno[,main.model,with=F],is.character)])
  
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


mgs.names.index <- grep("HG3A",names(pheno))
names(pheno)[mgs.names.index] <- cutlast(names(pheno)[mgs.names.index],9)

  # Exporting in STATA friendly format 


  variables.export <- c("SCAPISid", "ahi", "odi", "t90", "valid.ahi",
                      "visit.month",
                      main.model, dummy.names,
                      grep("HG3A",names(pheno),value=T))

  pheno <- pheno[valid.t90=="yes",]

  pheno[,leisurePA:=factor(leisurePA, labels=c("PA1","PA2","PA3","PA4"))]

  pheno[,educat:=factor(educat, labels=c("cat0", "cat1", "cat2", "cat3"))]


  require(foreign)
  write.dta(pheno[,variables.export,with=F], paste0(work,"pheno.dta"))
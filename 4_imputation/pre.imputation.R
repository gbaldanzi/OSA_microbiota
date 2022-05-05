# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script to prepare data to be used in STATA

# Because AHI have a lower sample size than ODI and T90, we are imputing 
# missing AHI values 

# STATA will be used to conducted the multiple imputation and the analyses with the 
# imputed data

# This script will prepare the data to be used at STATA

  rm(list=ls())

  library(data.table)
  library(vegan)


  # Folders 
  work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'

  # Import data
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))

  # Create dummy variables for factor variables 

  #Covariates 
  main.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","BMI")
  extended.model <- c(main.model,"Fibrer","Energi_kcal" ,"leisurePA", 
                      "educat","placebirth","month")
  
  # Complete cases for the extended model on the analysis with T90/ODI
  pheno <- pheno[valid.t90=="yes",]
  setnames(pheno,"visit.month","month")
  
  cc <- complete.cases(pheno[,extended.model,with=F])
  pheno <- pheno[cc,]
  

  
  # Rename levels for better handling at STATA
  pheno[,leisurePA := factor(leisurePA, levels(leisurePA), labels =c("PA1","PA2","PA3","PA4"))]
  pheno[,educat := factor(educat, levels(educat), labels = c("edu1","edu2","edu3", "edu4"))]
  pheno[month == "June.July", month := "Jun_Jul"]
  

  
  # Dummy variables
  temp.data <- pheno[,c("SCAPISid",extended.model),with=F]
  setDF(temp.data)
  rownames(temp.data) <- temp.data$SCAPISid
  temp.data <- as.data.frame(model.matrix(~.,temp.data[extended.model]))
  temp.data <- temp.data[,-which(names(temp.data)=="(Intercept)")]
  
  names.dummy.var <- names(temp.data)
  
  
  # Merge dummy variables 
  temp.data$SCAPISid <- rownames(temp.data)
  
  vars <- c("SCAPISid", "ahi", "t90", "odi","valid.ahi","valid.t90","shannon","WaistHip",grep("HG3A",names(pheno),value=T))
  
  pheno <- merge( pheno[, vars,with=F], temp.data, by="SCAPISid")
  
  pheno$monthJune.July <- NULL # all equal 0

# Makes species names shorter (some names were too long for STATA)

cutlast <- function(char,n){
  l <- nchar(char)
  a <- l-n+1
  return(substr(char,a,l))
}

  mgs.names.index <- grep("HG3A",names(pheno))
  names(pheno)[mgs.names.index] <- cutlast(names(pheno)[mgs.names.index],9)
  
  # Divide the data set in 4 equal sizes to run the imputation in parallel in STATA
  
  species.names <- grep("HG3A",names(pheno),value=T)
  
  species.data <- pheno[,c("SCAPISid",species.names),with=F]
  setDF(species.data)
  pheno0 <- pheno[,-which(names(pheno) %in% species.names),with=F]
  
  n <- round(length(species.names)/4)
  
  pheno_1 <- merge(pheno0, species.data[,c(1:n)],by="SCAPISid")
  pheno_2 <- merge(pheno0, species.data[,c(1,(n+1):(2*n))],by="SCAPISid")
  pheno_3 <- merge(pheno0, species.data[,c(1,(2*n+1):(3*n))],by="SCAPISid")
  pheno_4 <- merge(pheno0, species.data[,c(1,(3*n+1):ncol(species.data))],by="SCAPISid")

  # Exporting in STATA friendly format 
  require(foreign)
  
  write.dta(pheno_1, paste0(work,"pheno_1.dta"))
  write.dta(pheno_2, paste0(work,"pheno_2.dta"))
  write.dta(pheno_3, paste0(work,"pheno_3.dta"))
  write.dta(pheno_4, paste0(work,"pheno_4.dta"))
  
  
  
  
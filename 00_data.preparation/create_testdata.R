# Project Sleep Apnea and Gut Microbiota

# Gabriel Baldanzi 
# Create test data from a random sample of the real data and 
# applying random coefficient to modify the real data

# Not seed should be set to ensure that the test data cannot be linked to the real data

library(data.table)

rm(list=ls())

  # Folder 
  input <- "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/data_processed/"


  # Import data 
  pheno <-readRDS(paste0(input,"pheno_sleep_mgs.rds"))
  
  n <- 600

  pheno <- pheno[sample(1:n),] # random 600 participants 

  pheno$SCAPISid <- paste0("sample_",1:n) # substitute ids 
  
  
  # Antropometric data collection 
  
  pheno$AnthropometryCollectionDate <- pheno$AnthropometryCollectionDate + runif(n,-180,180) # round(1)
  pheno$AnthropometryCollectionDate <- pheno$AnthropometryCollectionDate + runif(n,-180,180) # round(2)
  

  # Continuous variables 

  cont.variables <- c("ahi","odi","t90","age","Alkohol","BMI","SBP_Mean","DBP_Mean",
                      "Fibrer","Energi_kcal","Hba1cFormattedResult", "HBFormattedResult",
                      "WaistHip")
  
  pheno.cont <- pheno[,cont.variables, with=F]
  
  pheno.cont <- apply(pheno.cont,2,function(x) x*runif(1,.8,1.2)) # multiple by a random coeff.
  
  colnames(pheno.cont) <- cont.variables
  
  
  # Categorical/Factor variables 
  # For simplicity, we reduced the number of levels in certain factor variables 
  factor.var <- c("Sex", "leisurePA","smokestatus", "o2utv4h", "BothFlO2utv4h",
                  "metformin", "hypermed", "dyslipmed",
                  "ppi", "lungdisease" , "atb6m",
                  "cqhe061", "cpap", "educat")
  
  pheno.factor <- apply(pheno[,factor.var,with=F],2,sample)
  
  colnames(pheno.factor) <- factor.var
  
  placebirth <-  sample(c("country1", "country2", "country3"), n, replace = T)
  plate <-  sample(c("plate1", "plate2", "plate3"), n, replace = T)
  

  
  
  # Species variables (30 species)
  species.names <- grep("HG3A",names(pheno),value=T)
  species.names <- sample(species.names,30) # Select 30 random species 
  species.matrix <- as.matrix(pheno[,species.names,with=F])
  
  random.coef1 <- runif(30,.8,1.2) 
  species.matrix <- t(random.coef1*t(species.matrix)) # Multiple by a random coefficient (round 1)
  
  random.coef2 <- runif(30,.8,1.2) 
  species.matrix <- t(random.coef2*t(species.matrix)) # Multiple by a random coefficient (round 2)
  
  species.t90 <- sample(grep("HG3A",names(pheno),value=T),5)
  species.t90 <- as.matrix(pheno[,species.t90,with=F])
  random.coef <- runif(5,.8,1.2) 
  species.t90 <- t(random.coef*t(species.t90))*(pheno$t90*.05)
  
  species.odi <- sample(grep("HG3A",names(pheno),value=T),5)
  species.odi <- as.matrix(pheno[,species.odi,with=F])
  random.coef <- runif(5,.8,1.2) 
  species.odi <- t(random.coef*t(species.odi))*(pheno$odi*.05)
  
  species.odit90 <- sample(grep("HG3A",names(pheno),value=T),10)
  species.odit90 <- as.matrix(pheno[,species.odit90,with=F])
  random.coef <- runif(10,.8,1.2) 
  species.odit90 <- t(random.coef*t(species.odit90))*(pheno$odi*.1)*(pheno$t90*.1)
  
  species.matrix <- cbind(species.matrixs, species.t90,species.odi,species.odit90)
  
  species.matrix <- t(apply(species.matrix,1,function(x) x/sum(x)))
  
  species.matrix <- species.matrix[sample(nrow(species.matrix)),]
  colnames(species.matrix) <- paste0("species____HG3A.",1:ncol(species.matrix))
  
  
  

  # Final data 
  pheno2 <- data.table( pheno$SCAPISid,pheno.factor,pheno$AnthropometryCollectionDate,
                       pheno.cont,species.matrix)

  saveRDS(pheno2, "pheno_sleep_mgs.rds")




                 
#simulate functional modules
  species.df <- data.table(colnames(species.matrix))
  species.df$module <- sample(c(1,1,2,2,2,3,3,3,3),ncol(species.matrix),replace=T) 

module1= species.df[]
module2= sample(names(mgs)[10:50],20)
module3= sample(names(mgs)[10:30],15)
module4= sample(colnames(species.matrix)[35:40],5)

module <- list(module1 = module1,
               module2 = module2,
               module3 = module3,
               module4=module4)


save(module,file=paste(output.folder,"/modules.RData",sep=""))
save(subpathway,file=paste(output.folder,"/metabolic_subpathways.RData",sep=""))

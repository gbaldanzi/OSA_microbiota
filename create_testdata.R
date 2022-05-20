# Project Sleep Apnea and Gut Microbiota

# Gabriel Baldanzi 
# Create test data from a random sample of the real data and 
# applying random coefficient to modify the real data

# Not seed should be set to ensure that the test data cannot be linked to the real data

library(data.table)
library(vegan)

rm(list=ls())

  # Folder 
  input <- "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/data_processed/"


  # Import data 
  pheno <-readRDS(paste0(input,"pheno_sleep_mgs.rds"))
  
  # Remove rare species
  species.names_all <- grep("HG3A",names(pheno),value=T)
  # presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
  data_pa <- decostand(x = pheno[,species.names_all, with=F], "pa")
  # calculate sum per species
  data_sum <- data.frame(prevalence = apply(data_pa, 2, sum)/nrow(data_pa))
  data_sum$MGS = rownames(data_sum)
  
  low.prevalence.species = data_sum$MGS[data_sum$prevalence <= 2/100] # 611 species have a prevalence > 2%
  pheno <- pheno[ , -low.prevalence.species, with=F] 
  
  
  # Reduce data to 600 participants 
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
  
  pheno.cont <- as.data.frame(pheno.cont)
  
  pheno.cont$AnthropometryCollectionDate <- pheno$AnthropometryCollectionDate
  
  
  # Categorical/Factor variables 
  # For simplicity, we reduced the number of levels in certain factor variables 
  factor.var <- c("Sex", "smokestatus",
                  "metformin", "hypermed", "dyslipmed", "diabmed",
                  "ppi", "lungdisease" , "atb6m",
                  "cqhe061", "cpap")
  
  pheno.factor <- apply(pheno[,factor.var,with=F],2,sample) # scramble factor variables 
  
  colnames(pheno.factor) <- factor.var
  
  pheno.factor <- as.data.frame(pheno.factor)
  
  pheno.factor$o2utv4h <- sample(pheno$o2utv4h)
  pheno.factor$BothFlO2utv4h <- sample(pheno$BothFlO2utv4h)
  pheno.factor$leisurePA <- sample(pheno$leisurePA)
  pheno.factor$educat <- sample(pheno$educat)
  pheno.factor$SCAPISid <- pheno$SCAPISid
  
  
  pheno.factor$placebirth <-  sample(c("country1", "country2", "country3"), n, replace = T)
  pheno.factor$plate <-  sample(c("plate1", "plate2", "plate3"), n, replace = T)
  

  
  
  # Species variables (30 species) - not associated with T90 or ODI
  species.names <- grep("HG3A",names(pheno),value=T)
  species.names <- sample(species.names,50) # Select 30 random species 
  species.matrix <- as.matrix(pheno[,species.names,with=F])
  
  random.coef1 <- rep(runif(30,.8,1.2),each=50)
  species.matrix <- random.coef1*species.matrix # Multiple by a random coefficient (round 1)
  
  random.coef2 <- rep(runif(30,.8,1.2),each=50)
  species.matrix <- random.coef2*species.matrix # Multiple by a random coefficient (round 2)
  
  species.matrix <- apply(species.matrix, 2, function(x) x+runif(n,.2,.5))
  species.matrix <- t(apply(species.matrix,1,function(x) x/sum(x)))
  
  
  species.matrix <- species.matrix[sample(nrow(species.matrix)),]
  
  
  # species associated with T90
  species.t90 <- species.matrix[,31:40]
  species.t90 <- apply(species.t90, 2 , function(x)  x[order(x)][rank(pheno.cont$t90)])
  index <- sample(1:(10*600),500)
  species.t90[index] <- sample(species.t90[index])
  
    # species associated with ODI 
  species.odi <- species.matrix[,41:50]
  species.odi <- apply(species.odi, 2 , function(x)  x[order(x)][rank(pheno.cont$odi)])
  index <- sample(1:(10*600),500)
  species.odi[index] <- sample(species.odi[index])

  
  species.matrix <- cbind(species.matrix[,1:30], species.t90, species.odi)

  colnames(species.matrix)[1:9] <- paste0("species____HG3A.000",1:9)
  colnames(species.matrix)[10:ncol(species.matrix)] <- paste0("species____HG3A.00",10:ncol(species.matrix))
  
  
  apply(species.t90, 2, cor, y = pheno.cont$t90, use = "complete.obs")
  
  
  # modules abundance 
  
  mod.abund <- sample(grep("MF0", names(pheno), value = T), 4)
  
  mod.abund <- apply(pheno[,mod.abund,with=F],2, function(x) x*runif(1,0.8,1.2))
  colnames(mod.abund) <- paste0("module",1:4)
  

  # Final data 
  pheno2 <- data.table(pheno.factor,
                       pheno.cont,species.matrix, mod.abund)
  
  setcolorder(pheno2, c("SCAPISid","Sex","age"))
  
  # work folder to save the final data
  if(!dir.exists("input")){dir.create("input")}
  
  output <- '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/sleepapnea_gut/'

  saveRDS(pheno2, paste0(output,"input/pheno_sleep_mgs.rds"))


  cutlast <- function(char,n){
    l <- nchar(char)
    a <- l-n+1
    return(substr(char,a,l))
  }


                 
#simulate functional modules
  species.df <- data.table(mgs = colnames(species.matrix))
  species.df$module <- sample(c(1,1,2,2,2,3,3,3,3),ncol(species.matrix),replace=T) 

module1= cutlast(species.df[module == 1 , mgs ],9)
module2= cutlast(species.df[module == 2 , mgs ],9)
module3= cutlast(species.df[module == 3 , mgs ],9)
module4= cutlast(sample(colnames(species.matrix)[35:50]),9)

MGS_HG3A.GMMs2MGS <- list(module1 = module1,
               module2 = module2,
               module3 = module3,
               module4 = module4)

output <- '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/sleepapnea_gut/'

save(MGS_HG3A.GMMs2MGS,file = paste0(output,"input/MGS_HG3A.GMMs2MGS.RData"))

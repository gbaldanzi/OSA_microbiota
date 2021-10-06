# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-09-21

# Last update: 2021-09-21

# Script to run the correlation between AHI and PC1-PC3
  
  library(data.table)
  library(fastDummies)
  library(ppcor)
  library(dplyr)
  library(tidyr)

  output = "/home/baldanzi/Sleep_apnea/Results/"
  input1 = '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/'

  # Import data and function 
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  mgs.pca_valid.ahi <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/mgs.pca_valid.ahi.rds")
  mgs.pca_valid.t90 <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/mgs.pca_valid.t90.rds")
  
  source(paste0(input1,"Spearman.correlation.function.R"))
  

  # Build function 
  
  cor_apnea_pc.fun <-  function(exposure, pc_data, pheno.data, covariates){
    
    dades <- copy(pheno.data)
    a <- c("SCAPISid",exposure,covariates)
    dades <- dades[,a,with=F]
    
    dades <- cbind(dades, pc_data$x[,c("PC1","PC2","PC3")])
  
    a= c("Sex")
    dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
    
    res <- spearman.function(x1=c("PC1","PC2","PC3"),x2=exposure,covari = covariates,data = dades)
    
    res = res[order(res$q.value),]
    
    names(res) <- c("outcome", "exposure", "rho", "p.value", "N", "method", "covariates","q.value")
    
    return(res)
    
  }
  
  clean.res.fun <- function(res,list.models){
   
   model.names <- names(list.models)
   
   for(i in 1:length(res)){  res[[i]]$model <- model.names[i]   }
  
   res.df <- do.call(rbind,res)
   rownames(res.df) <- NULL
   return(res.df)
  
  }
    
  # Models   
  model1 <- c("age", "Sex", "Alkohol","smokestatus")
  model2 <- c(model1,"BMI")
  model3 <- c(model2,"metformin","hypermed","dyslipmed","ppi","Fibrer",
              "Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")
    
  models <- list(model1 = model1, 
                      model2 = model2, 
                      model3 = model3)
  
  # Run correlations 
  
  res.pc.ahi <- lapply(models,cor_apnea_pc.fun,exposure="ahi",
                        pc_data = mgs.pca_valid.ahi,
                        pheno.data = pheno[valid.ahi=='yes',])
    
  
  res.pc.ahi <- clean.res.fun(res.pc.ahi, models)
    
  res.pc.t90 <- lapply(models,cor_apnea_pc.fun,exposure="t90",
                                 pc_data = mgs.pca_valid.t90,
                                 pheno.data = pheno[valid.t90=='yes',])
  
  res.pc.t90 <- clean.res.fun(res.pc.t90, models)
  
  res <- rbind(res.pc.ahi, res.pc.t90)
  
  fwrite(res, file=paste0(output,"cor_pc.tsv"), sep = '\t')
    
    
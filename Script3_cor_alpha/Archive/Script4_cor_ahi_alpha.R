# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: March, 2021
# Update: Sep 30, 2021

# This code will investigate the alpha diversity in relation AHI, OSAcat, and T90%.  
# Four models  and a sensitivity analysis will be used 
# Results will be saved in three files depending on the exposure (AHI,OSAcat,T90%)


# Correlation between AHI and Shannon ####
  
  # Prepare data 
  dades = copy(pheno[valid.ahi =="yes",])
  
  # Transforming two-level factor variables into numeric variables 
  a= c("Sex","hypermed","dyslipmed","ppi", "metformin")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
  
# Preparing exposure and  outcomes
  exposure="ahi"
  outcomes="shannon"

# Run Spearman correlation for the models.
  res.alpha = t(sapply(listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)}))

  res.alpha <-  as.data.frame(res.alpha)
  res.alpha$model <-  names(listmodels) 
  
  
  names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                       "N", "method", "covariates","model")
  
  fwrite(res.alpha, file = paste0(output,"cor_ahi_alpha.tsv"), sep="\t")
  
  res.alpha.ahi <- res.alpha


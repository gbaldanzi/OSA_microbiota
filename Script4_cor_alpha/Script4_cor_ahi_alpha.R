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
  res.alpha$model <-  c("model1", "model2", "model3","SA2")


#----------------------------------------------------------------------------#
# Sensitivity analysis 1
  
  # Prepare data 
  dades.sa = copy(pheno[valid.ahi =="yes",])
  dades.sa <- dades.sa[ppi=="no",]
  dades.sa <- dades.sa[metformin=="no",]
  dades.sa <- dades.sa[hypermed=="no",]
  dades.sa <- dades.sa[dyslipmed=="no",]
  
  # Transforming two-level factor variables into numeric variables 
  a= c("Sex")
  dades.sa[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades.sa[,a, with=F]))))]

  # Spearman correlation 
  res.sa = spearman.function(x1=outcomes,x2=exposure,covari = SA,data = dades.sa)
  res.sa$model= "SA"

  # Naming columns
  res.alpha = rbind(res.alpha,res.sa)
  
  names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","model")

  fwrite(res.alpha, file = paste0(output,"cor_ahi_alpha.tsv"), sep="\t")

#-----------------------------------------------------------------------------#

# Correlation between AHI and Shannon by BMI group ####

# the correlation between AHI and shannon will be separately by the 3 BMI groups:
# <25, >=25 & <30, and >=30.
  
  setDT(dades)
  setDT(dades.sa)

  # By BMI group 
  for(group in unique(pheno[,BMIcat])){
    dades2 = dades[BMIcat==group,]


  # Run Spearman correlation for the models.
  res.alpha = sapply(listmodels,function(mod){
     spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades2)})
  res.alpha <- as.data.frame(t(res.alpha))
  
  res.alpha$model = c("model1", "model2", "model3","SA2")
  
  res.alpha$bmi= group


  
#----------------------------------------------------------------------------#
  message("Sensitivity analysis")
  # Sensitivity analysis 

  # BMI group
  dades2 <- dades.sa[dades.sa$BMIcat==group,]

  # Spearman correlation 
  res.sa = spearman.function(x1=outcomes,x2=exposure,covari = SA,data = dades2)

  res.sa$model= "SA"
  res.sa$bmi= group

  res.alpha = rbind(res.alpha,res.sa)
  
  #naming columns
  names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","model","bmi")

  fwrite(res.alpha, file = paste0(output,"cor_ahi_alpha_bmi",group,".tsv"), sep="\t")
  }


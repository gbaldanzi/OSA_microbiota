# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-16

# This code will investigate the association between species and AHI 



# Correlation between AHI and MGS - Step1 ####

#Preparing exposure, outcomes, and covariates
exposure="ahi"
outcomes=grep("___",names(dades),value=T)

# Running correlation 
  res<-   spearman.function(x1=outcomes,x2=exposure,covari = basic.model,data = dades)

  res.ahi <- res


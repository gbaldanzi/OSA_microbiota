# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-09-16

# Inferential Statistics 

# This code will investigate the association between MGS and AHI 

# Sensitivity Analysis 3 - additional adjustment for BMI

# Results saved at "/home/baldanzi/Sleep_apnea/Results/corSA3_ahi_mgs.tsv"


  # Importing data
  dades <- pheno[valid.ahi=="yes",] 

  # Correlation between AHI and MGS - Sensitivity analysis 3 ####

  #Preparing exposure, outcomes, and covariates
  exposure="ahi"

  # Running correlation 
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = SA3,data = dades)

  res = res[order(res$q.value),]

  res$model= "SA3"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")

  fwrite(res, file = paste0(output,"corSA3_ahi_mgs.tsv"), sep="\t")
  

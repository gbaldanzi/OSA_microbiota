# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-01

# Inferential Statistics 

# This code will investigate the association between MGS and AHI 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/corsa_ahi_mgs.tsv"

# Importing data
dades <- copy(pheno[valid.ahi=='yes',])

# Correlation between AHI and MGS  ####

#Preparing exposure, outcomes, and covariates
exposure="ahi"

# Running correlation 
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = model2,data = dades)

  res = res[order(res$q.value),]

  res$model= "SA"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")

  fwrite(res, file = paste0(output,"corsa_ahi_mgs.tsv"), sep="\t")

  res.ahi <- res


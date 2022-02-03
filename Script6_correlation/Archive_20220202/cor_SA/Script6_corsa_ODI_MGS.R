# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-12-15


# This code will investigate the association between MGS and ODI in a sensitivity analysis
#excluding participants that have used PPI, metformin, self-reported anti-hypertensive or 
#cholesterol lowering medication

# Results saved at "/home/baldanzi/Sleep_apnea/Results/corsa_odi_mgs.tsv"


# Importing data
  dades <- copy(pheno[valid.t90=='yes',])


# Correlation between t90 and MGS  ####

#Preparing exposure, outcomes, and covariates
exposure="odi"


# Running correlation 
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = model2,data = dades)

  res = res[order(res$q.value),]

  res$model= "SA"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")

fwrite(res, file = paste0(output,"corsa_odi_mgs.tsv"), sep="\t")

  res.odi <- res

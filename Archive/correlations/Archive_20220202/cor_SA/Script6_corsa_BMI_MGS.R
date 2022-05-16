# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-01

# Inferential Statistics 

# This code will investigate the association between MGS and BMI 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/corsa_BMI_mgs.tsv"

# Transforming two-level factor variables into numeric variables 
dades = copy(pheno)


# Correlation between BMI and MGS  ####

#Preparing exposure
exposure="BMI"


# Running correlation 
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = model2,data = dades)

  res = res[order(res$q.value),]

  res$model= "SA"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")

fwrite(res, file = paste0(output,"corsa_BMI_mgs.tsv"), sep="\t")

 res.bmi <- res
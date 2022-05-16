# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-16

# Inferential Statistics 

# This code will investigate the association between BMI and MGS 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/cor_BMI_mgs.tsv"

# Transforming two-level factor variables into numeric variables 
dades = copy(pheno)
a= c("Sex", "metformin","hypermed","dyslipmed","ppi")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]


# Correlation between BMI and MGS - Basic model ####

#Preparing exposure, outcomes, and covariates
exposure="BMI"
outcomes=grep("___",names(pheno),value=T)

#Covariates 
basic.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")

# Running correlation 
  message("Correlation between MGS and BMI - Step1")
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = basic.model,data = dades)

  res = res[order(res$q.value),]

  res.bmi <- res


# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Update: Dec 14, 2021

# This code will investigate the alpha diversity in relation ODI.  


# This code must be preceded by Script4_cor_main.R

# Correlation between ODI and Shannon ####
  
  # Prepare data 
  dades = copy(pheno[valid.t90=='yes',])
  
  # Transforming two-level factor variables into numeric variables 
  a= c("Sex","hypermed","dyslipmed","ppi", "metformin")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
  
# Preparing exposure and  outcomes
  exposure="odi"
  outcomes="shannon"

# Run Spearman correlation for the models.
  res.alpha = t(sapply(listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)}))

  res.alpha <-  as.data.frame(res.alpha)
  res.alpha$model <-  names(listmodels) 


  names(res.alpha) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","model")

  fwrite(res.alpha, file = paste0(output,"cor_odi_alpha.tsv"), sep="\t")

  res.alpha.odi <- res.alpha

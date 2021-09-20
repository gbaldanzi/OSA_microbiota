# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-17

# Inferential Statistics 

# This code will investigate the association between MGS and BMI 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/cor2_BMI_mgs.tsv"

# MGSs identified in model 1
#mgs.m1 <- readRDS('/home/baldanzi/Sleep_apnea/Results/mgs.m1.rds')

# Transforming two-level factor variables into numeric variables 
dades = copy(pheno)
a= c("Sex", "metformin","hypermed","dyslipmed","ppi")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Correlation between BMI and MGS - Step2 ####

  #Preparing exposure, outcomes, and covariates
  exposure="BMI"
  outcomes <- mgs.m1

  # Running correlation 
  message("Correlation between MGS and BMI - Step2")
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = model2,data = dades)

  res = res[order(res$q.value),]

  res$model= "model2"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
               "N", "method", "covariates","q.value","model")

  fwrite(res, file = paste0(output,"cor2_BMI_mgs.tsv"), sep="\t")
  #fwrite(res, file = paste0(output,"cor2_BMI_mgs_filter001.tsv"), sep="\t")
  
  res.bmi <- res

#--------------------------------------------------------------------------#
# Merging results with taxonomy information #### 

#taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
#setnames(taxonomy,"maintax_mgs","MGS")

#dades <- fread(paste0(output,"cor2_BMI_mgs.tsv"))
#dades <- merge(dades, taxonomy, by="MGS", all.x=T)
#fwrite(dades, file=paste0(output,"cor2_BMI_mgs.tsv"))
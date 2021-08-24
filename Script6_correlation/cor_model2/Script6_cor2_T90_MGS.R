# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-17

# Inferential Statistics 

# This code will investigate the association between MGS and t90 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/cor2_t90_mgs.tsv"

# Transforming two-level factor variables into numeric variables 
  dades = copy(pheno[valid.t90=="yes",])
  a= c("Sex", "metformin","hypermed","dyslipmed","ppi")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]


# Correlation between t90 and MGS - Step2 ####

#Preparing exposure, outcomes, and covariates
  exposure <- "t90"
  outcomes <-  mgs.m1

# Running correlation 
  message("Correlation between MGS and T90 - Step2")
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = model2,data = dades)

  res = res[order(res$q.value),]

  res$model= "model2"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
               "N", "method", "covariates","q.value","model")

#fwrite(res, file = paste0(output,"cor2_t90_mgs.tsv"), sep="\t")
  fwrite(res, file = paste0(output,"cor2_t90_mgs_filter001.tsv"), sep="\t")

#--------------------------------------------------------------------------#
# Merging results with taxonomy information #### 

#taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
#setnames(taxonomy,"maintax_mgs","MGS")

#dades <- fread(paste0(output,"cor2_t90_mgs.tsv"))
#dades <- merge(dades, taxonomy, by="MGS", all.x=T)
#fwrite(dades, file=paste0(output,"cor2_t90_mgs.tsv"))
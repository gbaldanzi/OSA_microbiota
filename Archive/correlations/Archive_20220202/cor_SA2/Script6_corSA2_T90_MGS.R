# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-10

# Inferential Statistics 

# This code will investigate the association between MGS and t90 

# Sensitivity Analysis 2 - additional adjustment for sleep time 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/corSA2_t90_mgs.tsv"

  # Prepare data
  dades <- pheno[valid.t90=="yes",] 

  # Correlation between t90 and MGS - Sensitivity analysis 2 ####

  #Preparing exposure, outcomes, and covariates
  exposure="t90"
  

  # Running correlation 
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = SA2,data = dades)

  res = res[order(res$q.value),]

  res$model= "SA2"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
               "N", "method", "covariates","q.value","model")

  fwrite(res, file = paste0(output,"corSA2_t90_mgs.tsv"), sep="\t")

#--------------------------------------------------------------------------#
# Merging results with taxonomy information #### 

taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
setnames(taxonomy,"maintax_mgs","MGS")

dades <- fread(paste0(output,"corSA2_t90_mgs.tsv"))
dades <- merge(dades, taxonomy, by="MGS", all.x=T)
fwrite(dades, file=paste0(output,"corSA2_t90_mgs.tsv"))
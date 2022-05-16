# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-10

# Inferential Statistics 

# This code will investigate the association between MGS and AHI 

# Sensitivity analysis 2 - additional adjustment for sleeptime 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/corSA2_ahi_mgs.tsv"


  # Importing data
  dades <- pheno[valid.ahi=="yes",] 

  # Correlation between AHI and MGS - Sensitivity analysis 2 ####

  #Preparing exposure, outcomes, and covariates
  exposure="ahi"

  # Running correlation 
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = SA2,data = dades)

  res = res[order(res$q.value),]

  res$model= "SA2"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")

  fwrite(res, file = paste0(output,"corSA2_ahi_mgs.tsv"), sep="\t")
  
  #--------------------------------------------------------------------------#
# Merging results with taxonomy information #### 

taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
setnames(taxonomy,"maintax_mgs","MGS")

dades <- fread(paste0(output,"corSA2_ahi_mgs.tsv"))
dades <- merge(dades, taxonomy, by="MGS", all.x=T)
fwrite(dades, file=paste0(output,"corSA2_ahi_mgs.tsv"))
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-17

# Inferential Statistics 

# This code will investigate the association between MGS and t90 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/cor2_t90_mgs.tsv"

# input and output folder 
input = "/home/baldanzi/Sleep_apnea/Results/"

# MGSs identified in model 1
#mgs.m1 <- readRDS('/home/baldanzi/Sleep_apnea/Results/mgs.m1.rds')

# Importing data
  valid.t90 <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/valid.t90_MGS.shannon_Upp.rds")


# Transforming two-level factor variables into numeric variables 
  dades = copy(valid.t90)
  a= c("Sex", "metformin","hypermed","dyslipmed","ppi")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming characters to factor variables 
  dades[,plate:=as.factor(dades$plate)]

# Correlation between t90 and MGS - Step1 ####

#Preparing exposure, outcomes, and covariates
  exposure="t90"
  m = grep("____",names(dades),value=T)
  m = m[m %in% mgs.m1]
  outcomes=m

# Running correlation 
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
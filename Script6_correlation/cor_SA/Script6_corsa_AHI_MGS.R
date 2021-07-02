# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-01

# Inferential Statistics 

# This code will investigate the association between MGS and AHI 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/corsa_ahi_mgs.tsv"

# input and output folder 
input = "/home/baldanzi/Sleep_apnea/Results/"

# MGSs identified in model 1
mgs.m1 <- readRDS('/home/baldanzi/Sleep_apnea/Results/mgs.m1.rds')

# Importing data
valid.ahi <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/validsleep_MGS.shannon_Upp.rds")

#Calculating  MGS prevalence ####
noms=grep("____",names(valid.ahi),value=T)
# presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
data_pa <- decostand(x = valid.ahi[,noms,with=F], "pa")

# calculate sum per species
data_sum <- data.frame(prevalence=apply(data_pa, 2, sum))
data_sum$MGS = rownames(data_sum)
a = data_sum$MGS[data_sum$prevalence<5] #19 MGS are only present in less than 5 individuals

# Removing MGS that are rare
valid.ahi <- valid.ahi[ , -a, with=F] 

# Remove medication users
valid.ahi <-  valid.ahi[ppi == "no",] #60
valid.ahi <-  valid.ahi[metformin == "no",] #57 
valid.ahi <-  valid.ahi[hypermed == "no",] #593
valid.ahi <-  valid.ahi[dyslipmed == "no",] # 239

# Transforming two-level factor variables into numeric variables 
dades = copy(valid.ahi)
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming characters to factor variables 
dades[,plate:=as.factor(dades$plate)]


# Correlation between AHI and MGS  ####

#Preparing exposure, outcomes, and covariates
exposure="ahi"
m = grep("____",names(dades),value=T)
m = m[m %in% mgs.m1]
outcomes=m

# Running correlation 
  res <-   spearman.function(x1=outcomes,x2=exposure,covari = model2,data = dades)

  res = res[order(res$q.value),]

  res$model= "SA"

  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")

fwrite(res, file = paste0(output,"corsa_ahi_mgs.tsv"), sep="\t")

#--------------------------------------------------------------------------#
# Merging results with taxonomy information #### 

taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
setnames(taxonomy,"maintax_mgs","MGS")

dades <- fread(paste0(output,"corsa_ahi_mgs.tsv"))
dades <- merge(dades, taxonomy, by="MGS", all.x=T)
fwrite(dades, file=paste0(output,"corsa_ahi_mgs.tsv"))
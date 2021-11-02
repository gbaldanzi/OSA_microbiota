# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-24

model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
model2 <- c(model1,"metformin","hypermed","dyslipmed","ppi","Fibrer",
            "Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")

# Loading packages 
pacman::p_load(data.table, ppcor, fastDummies, vegan)

# Output folders 
input1 = '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/'
output = "/home/baldanzi/Sleep_apnea/Results/"

# Importing data
pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

# Removing rare MGS   
  # Calculating MGS prevalence 
  noms=grep("____",names(pheno),value=T)
# presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
  data_pa <- decostand(x = pheno[,noms,with=F], "pa")
# calculate sum per species
  data_sum <- data.frame(prevalence=apply(data_pa, 2, sum),
                       percentage = apply(data_pa,2,sum)/nrow(pheno))
  data_sum$MGS = rownames(data_sum)
  a = data_sum$MGS[data_sum$prevalence<5] #19 MGS are only present in less than 5 individuals
  #a = data_sum$MGS[data_sum$percentage<.01] #383 MGS are present in less than 1% of individuals

  pheno <- pheno[ , -a, with=F] 
  
  
# Correlations 
    source(paste0(input1,"Spearman.correlation.function.R")) # Correlation function 
  
    source('cor_model1/Script6_cor1_AHI_MGS.R')
    source('cor_model1/Script6_cor1_T90_MGS.R')
    source('cor_model1/Script6_cor1_BMI_MGS.R')
  
  res <- rbind(res.ahi, res.t90, res.bmi)
  
  res$model= "model1"
  
  names(res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                 "N", "method", "covariates","q.value","model")
  
  # Merging results with taxonomy information #### 
  
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  setnames(taxonomy,"maintax_mgs","MGS")
  
  res <- merge(res, taxonomy, by="MGS", all.x=T)
  
  fwrite(res, file = paste0(output,"cor_all.var_mgs.tsv"))
  
  res$q.value[res$q.value>=0.001] <- round(res$q.value[res$q.value>=0.001] , digits = 3)
  mgs.m1 <- unique(res[res$q.value<.05,"MGS"]) 
    
  saveRDS(mgs.m1, paste0(output,'mgs.m1.rds'))
  

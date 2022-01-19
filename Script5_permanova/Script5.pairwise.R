# Script 5 - Main 

# Script to run the permutation analysis between OSA/AHI/T90 and BC/AD

# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)

# Import data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

# Transforming two-level factor variables into numeric variables 
  a= c("Sex","ppi","metformin","hypermed","dyslipmed")
  pheno[,(a):=as.data.frame(data.matrix(data.frame(unclass(pheno[,a, with=F]))))]

#Covariates 
  # basic.model : adjust for age + sex + alcohol + smoking + plate + received 
  basic.model<-   c("age", "Sex", "Alkohol","smokestatus","plate")

  # full.model = basic.model + fiber intake + Energy intake + physical activity + education + country of birth + ppi + metformin +  antihypertensive + cholesterol-lowering 
  full.model <-  c(basic.model, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", "visit.month", "metformin","hypermed","dyslipmed","ppi")
  

  
  source('Script5_permanova/perma_ahi_bc/perma_pairwise_main.R')

  source('Script5_permanova/perma_odi_bc/perma_pairwise_main.R')

  source('Script5_permanova/perma_t90_bc/perma_pairwise_main.R')
  
  
 

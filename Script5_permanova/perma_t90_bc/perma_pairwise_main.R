# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: 2021-12-22
# Last Update: - 2021-12-22

# This code will run pairwise comparisons between T90 severity group in 
# relation to the beta-diversity (Bray Curtis Dissimilarity) 

# Analysis are run using a PERMANOVA approach


# Loading packages 
  pacman::p_load(data.table, vegan,parallel)


  output = "/home/baldanzi/Sleep_apnea/Results/"
  

  # Importing data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

  # Transforming two-level factor variables into numeric variables 
  a= c("Sex","ppi","metformin","hypermed","dyslipmed")
  pheno[,(a):=as.data.frame(data.matrix(data.frame(unclass(pheno[,a, with=F]))))]

  # Importing BC matrix 
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/T90.BCmatrix.csv', header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)

  source('perma.pairwise.fun.R')


# Making sure that BC and dataset have the same order of observations 
  pheno[,valid.t90="yes"]
  pheno <-  pheno[match(rownames(BC),pheno$SCAPISid),]
  

#Covariates 
  # model 1 : adjust for age + sex + alcohol + smoking + plate + received 
  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate")
  # model 2 = model 1 + BMI 
  model2 <-  c(model1,"BMI")
  # model 3 = model 2 + fiber intake + Energy intake + physical activity + education + country of birth + ppi + metformin +  antihypertensive + cholesterol-lowering 
  model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth","visit.month","metformin","hypermed","dyslipmed","ppi")


# Runing PERMANOVA in parallel ####
  
  pheno[OSAcat=="no OSA", OSAcat:="no_OSA"]
  
  list.res = vector(mode = "list",length = 0)
  
  g1 = "no_OSA" ; g2 = c("Mild","Moderate","Severe")
  
  for(i in 1:3){
  message(paste(g1,g2[i],sep = "_"))
  list.res[[paste(g1,g2[i],sep = "_")]] <-  pairwise.perma.fun(outcome = "BC", group1 = g1, group2 = g2[i], covari = model3, data = pheno, nodes = 16)
  }
  
  g1 = "Mild" ; g2 = c("Moderate","Severe")
  
  for(i in 1:2){
    message(paste(g1,g2[i],sep = "_"))
  list.res[[paste(g1,g2[i],sep = "_")]] <-  pairwise.perma.fun(outcome = "BC", group1 = g1, group2 = g2[i], covari = model3, data = pheno, nodes = 16)
  }
  
  g1 = "Moderate" ; g2 = c("Severe")
  
  for(i in 1:1){
    message(paste(g1,g2[i],sep = "_"))
  list.res[[paste(g1,g2[i],sep = "_")]] <-  pairwise.perma.fun(outcome = "BC", group1 = g1, group2 = g2[i], covari = model3, data = pheno, nodes = 16)
  }
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results.rds"))
#---------------------------------------------------------------------------#

  # Produce a final summary results table with FDR-p-values 
  
  list.res = readRDS("/home/baldanzi/Sleep_apnea/Results/pairwise.perma.results.rds")
  
  list.res$table <- clean.res(list.res[1:6]) 
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results.rds"))
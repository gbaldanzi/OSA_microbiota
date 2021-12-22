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
  pheno <- pheno[valid.t90=="yes",]

  # Transforming two-level factor variables into numeric variables 
  a= c("Sex","ppi","metformin","hypermed","dyslipmed")
  pheno[,(a):=as.data.frame(data.matrix(data.frame(unclass(pheno[,a, with=F]))))]

  # Importing BC matrix 
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/T90.BCmatrix.csv', header=T, sep = ',')
  BC = as.matrix(BC)
  row.names(BC) = colnames(BC)

  source('Script0_functions/perma.pairwise.fun.R')


# Making sure that BC and dataset have the same order of observations 
  pheno <-  pheno[match(rownames(BC),pheno$SCAPISid),]
  

#Covariates 
  # model 1 : adjust for age + sex + alcohol + smoking + plate + received 
  model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate")
  # model 2 = model 1 + BMI 
  model2 <-  c(model1,"BMI")
  # model 3 = model 2 + fiber intake + Energy intake + physical activity + education + country of birth + ppi + metformin +  antihypertensive + cholesterol-lowering 
  model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth","visit.month","metformin","hypermed","dyslipmed","ppi")


# Runing PERMANOVA in parallel ####
  
  pheno[, t90cat:=factor(t90cat, levels(t90cat), 
                         labels = c("t0", "t1", "t2", "t3"))]
  
  
  list.res <- pairwise.perma.fun(outcome="BC", group_var="t90cat", covari=model3, data=pheno, nodes=16)
  
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results_t90.rds"))
#---------------------------------------------------------------------------#

  # Produce a final summary results table with FDR-p-values 
  
  list.res = readRDS("/home/baldanzi/Sleep_apnea/Results/pairwise.perma.results.rds")
  
  list.res$table <- clean.res(list.res[1:6]) 
  
  saveRDS(list.res,file=paste0(output,"pairwise.perma.results.rds"))
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlation between Alpha diversity (Shannon index) and the variables AHI, T90 and ODI 
# Update: 2022-02-02

# Defining output folders 
  input <-  "/home/baldanzi/Datasets/sleep_SCAPIS/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"

# Importing data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  

#Spearman correlation function ####
source('Script0_functions/Spearman.correlation.function.R')

# Models ####

  #Covariates 
  # basic.model : adjust for age + sex + alcohol + smoking + plate + received 
  basic.model <-   c("age", "Sex", "Alkohol","smokestatus","plate")
  # full model = basic model + fiber intake + Energy intake + physical activity + education + country of birth + ppi + metformin +  antihypertensive + cholesterol-lowering 
  full.model <-  c(basic.model, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", "visit.month", "metformin","hypermed","dyslipmed","ppi")
  # full.model_bmi = full model + BMI 
  full.model_bmi = c(full.model, "BMI")
  
# Run correlation analysis 
  
  models = list(basic.model, full.model, full.model_bmi)
  exposures = c("ahi","t90","odi")
  
  res.alpha = vector(mode = "list")
  
  for(exp in exposures){
    
    res.temp <- lapply(models,spearman.function,
                     x1="shannon",x2=exp,data = pheno)
    res.alpha[[exp]] <- do.call(rbind,res.temp)
  
  }

  res.alpha <- do.call(rbind,res.alpha)
  
# Save results 
  
  fwrite(res.alpha, file= paste0(output, "cor_all.var_alpha.tsv"))
  
  
  
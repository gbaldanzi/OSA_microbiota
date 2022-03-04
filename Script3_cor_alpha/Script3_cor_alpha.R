# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlation between Alpha diversity (Shannon index) and the variables AHI, T90 and ODI 
# Update: 2022-02-02

# Defining output folders 
  input <-  "/home/baldanzi/Datasets/sleep_SCAPIS/"
  results.folder <-  "/home/baldanzi/Sleep_apnea/Results/"

# Importing data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  

#Spearman correlation function ####
source('Script0_functions/Spearman.correlation.function.R')

# Models ####

  #Covariates 
  main.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","BMI")
  
  extended.model <- c(main.model,"Fibrer","Energi_kcal" ,"leisurePA", 
                  "educat","placebirth","visit.month")
   
  medication.model <- c(main.model, "metformin","hypermed", "dyslipmed", "ppi")

# Run correlation analysis 
  
  exposures = c("ahi","t90","odi")
  
  res.alpha = vector(mode = "list")
  
  
    res.temp1 <- lapply(exposures,spearman.function,
                     x1="shannon",covari=main.model,data = pheno)
    
    res.temp2 <- lapply(exposures,spearman.function,
                        x1="shannon",covari=extended.model,data = pheno)
    
    res.temp3 <- lapply(exposures,spearman.function,
                        x1="shannon",covari=medication.model,data = pheno)
    
    res.alpha <- do.call(rbind,c(res.temp1,res.temp2,res.temp3))
  
# Save results 
  
  fwrite(res.alpha, file= paste0(results.folder, "cor_all.var_alpha.tsv"))
  
  
  
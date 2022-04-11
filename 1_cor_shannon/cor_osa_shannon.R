# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Correlation between Alpha diversity (Shannon index) and the variables AHI, T90 and ODI 

library(data.table)

  # Defining folders 
  input <-  "/home/baldanzi/Datasets/sleep_SCAPIS/"
  results.folder <-  "/home/baldanzi/Sleep_apnea/Results/"

  # Importing data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  

  #Spearman correlation function ####
  source('0_functions/Spearman.correlation.function.R')

  # Models ####

  #Covariates 
  main.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","BMI")
  
  extended.model <- c(main.model,"Fibrer","Energi_kcal" ,"leisurePA", 
                  "educat","placebirth","visit.month")
   
  # Spearman correlation  
  
  exposures = c("ahi","t90","odi")
  
  res.alpha = vector(mode = "list")
  
  
    res.main.model <- lapply(exposures,spearman.function,
                     x1="shannon",covari=main.model,data = pheno)
    
    res.ext.model <- lapply(exposures,spearman.function,
                        x1="shannon",covari=extended.model,data = pheno)
    
    res.alpha <- do.call(rbind,c(res.main.model,res.ext.model))
  
# Save results 
  
  fwrite(res.alpha, file= paste0(results.folder, "cor_all.var_alpha.tsv"))
  
  
  
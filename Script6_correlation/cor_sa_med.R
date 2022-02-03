# Script 6 - Sensitivity analysis 

# Gabriel Baldanzi 
# 2021-09-02

# Last update: 2022-02-102

 message("Sensitivity Analysis 1")

# Sensitivity analysis 
# Remove users of ppi, metformin, hypertensive and cholesterol-lowering medication 
#  adjust for basic.model + place birth + education + leisure PA + fiber + total energy intake 
  
  # Import MGS identified in model 1
  # mgs.m1  = readRDS(paste0(output,'mgs.m1.rds'))
  
  # Importing data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  
  # Remove medication users
  pheno <-  pheno[ppi == "no",] #60
  pheno <-  pheno[metformin == "no",] #57 
  pheno <-  pheno[hypermed == "no",] #593
  pheno <-  pheno[dyslipmed == "no",] # 239
  
  basic.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  SA1 <- c(model1,"Fibrer","Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")
  
  outcomes <- mgs.m1


# Correlation scripts  
  
  outcomes=grep("___",names(dades),value=T)
  exposures <- c("ahi","t90","odi","BMI")
  
  res.sa1 <- lapply(exposures,spearman.function, 
                            x1=outcomes,
                            covari = SA1,
                            data = pheno)
  
  res.sa1 <- do.call(rbind,res.sa1)
  
  setDT(res.sa1)
  setnames(res.sa1,"x","MGS")
  

  fwrite(res.sa1, file = paste0(output,"corsa_all.var_mgs.tsv"))
  

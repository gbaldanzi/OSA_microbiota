# Script 6 - Sensitivity analysis 

# Gabriel Baldanzi 

# Last update: 2022-02-15

# Sensitivity analysis 
# Remove users of ppi, metformin, hypertensive and cholesterol-lowering medication 
#  adjust for basic.model + place birth + education + leisure PA + fiber + total energy intake 

  # Importing data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  
  # Remove medication users
  pheno <-  pheno[ppi == "no",] #60
  pheno <-  pheno[metformin == "no",] #57 
  pheno <-  pheno[hypermed == "no",] #593
  pheno <-  pheno[dyslipmed == "no",] # 239
  
  SAmed <- c(basic.model,"BMI","Fibrer","Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")
  
  # Importing species that were associated with the exposures in the full model 
  # mgs.m1 <- readRDS(paste0(results.folder, 'mgs.m1.rds'))
  
# Correlation scripts  
  
  # mgs.m1 <- readRDS(paste0(results.folder,'mgs.m1.rds'))
  outcomes <- mgs.m1
  exposures <- c("ahi","t90","odi")
  
  res.sa.med <- lapply(exposures,spearman.function, 
                            x1=outcomes,
                            covari = SAmed,
                            data = pheno)
  
  res.sa.med <- do.call(rbind,res.sa.med)
  
  setDT(res.sa.med)
  setnames(res.sa.med,"x","MGS")
  

  fwrite(res.sa.med, file = paste0(results.folder,"corsamed_all.var_mgs.tsv"))
  

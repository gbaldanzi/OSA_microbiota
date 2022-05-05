# Project sleep apnea and gut microbiota 

# Gabriel Baldanzi 

# Association between T90/ODI and species stratified by the hemoglobin

# In this script, we will only include the species that were FDR significant in the main model 

# This sensitivity analysis is adjusted for all main model covariates


  library(data.table)
  library(dplyr)
  library(tidyr)

  # Folders 
  results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
  work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'

  # Spearman with boostrap function
  source("0_functions/cor.boot.fun.R")
  
  # Import data
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))


  # Species associated with T90 or ODI in the extended
  res.main.model <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
  
  mgs.t90 <- res.main.model[q.value<.05 & exposure=="t90",MGS]
  mgs.odi <- res.main.model[q.value<.05 & exposure=="odi",MGS]
  
  
  # Model covariates 
  extended.model <- c("age", "Sex", "Alkohol","smokestatus",
                  "plate","BMI", "Fibrer","Energi_kcal" ,"leisurePA", 
                  "educat","placebirth","visit.month")
  
  
  # HB groups based on Sex-specific median Hb level 
  median.m <- median(pheno[Sex=="male",HBFormattedResult],na.rm=T)
  median.f <- median(pheno[Sex=="female",HBFormattedResult],na.rm=T)
  
  pheno[Sex == "male" & HBFormattedResult < median.m, Hbgroup := "low"]
  pheno[Sex == "female" & HBFormattedResult < median.f, Hbgroup := "low"]
  pheno[Sex == "male" & HBFormattedResult >= median.m, Hbgroup := "high"]
  pheno[Sex == "female" & HBFormattedResult >= median.f, Hbgroup := "high"]
  pheno[,Hbgroups := factor(Hbgroup, levels = c("low","high"))]
  
  
  # Restrict data to those with valid T90/ODI
  pheno <- pheno[valid.t90=="yes",]
  
  # T90 analysis
  ## Low HB
  message("T90 analysis - low HB")
  res.t90.hb.low <- lapply(mgs.t90, cor.boot, x = "t90", z=extended.model, data=pheno[Hbgroup =="low",])
  res.t90.hb.low <- do.call(rbind, res.t90.hb.low)
  res.t90.hb.low$Hbgroup <- "low"
  
  ## High HB
  message("T90 analysis - high HB")
  res.t90.hb.high <- lapply(mgs.t90, cor.boot, x = "t90", z=extended.model, data=pheno[Hbgroup =="high",])
  res.t90.hb.high <- do.call(rbind, res.t90.hb.high)
  res.t90.hb.high$Hbgroup <- "high"
  
  res.t90 <- rbind(res.t90.hb.high, res.t90.hb.low)
  
  res.t90 <- pivot_wider(res.t90, id_cols = c(outcome,exposure), names_from = Hbgroup, 
                           values_from = c(rho,se,conf.int, p.value,N, covariates))
  
  
  # ODI analysis 
  ## Low HB
  message("ODI analysis - low HB")
  res.odi.hb.low <- lapply(mgs.odi, cor.boot, x = "odi", z=extended.model, data=pheno[Hbgroup =="low",])
  res.odi.hb.low <- do.call(rbind, res.odi.hb.low)
  res.odi.hb.low$Hbgroup <- "low"
  
  ## High HB
  message("ODI analysis - high HB")
  res.odi.hb.high <- lapply(mgs.odi, cor.boot, x = "odi", z=extended.model, data=pheno[Hbgroup =="high",])
  res.odi.hb.high <- do.call(rbind, res.odi.hb.high)
  res.odi.hb.high$Hbgroup <- "high"
  
  res.odi <- rbind(res.odi.hb.high, res.odi.hb.low)
  
  res.odi <- pivot_wider(res.odi, id_cols = c(outcome,exposure), names_from = Hbgroup, 
                         values_from = c(rho,se,conf.int, p.value,N, covariates))
  
  # Heterogeneity test ####
  
  heterog.test.fun <- function(res){
  
  diff <-   res$rho_low - res$rho_high
  se <- sqrt((res$se_low^2)+(res$se_high^2))
  z.score <- abs(diff/se)
  pv <- pnorm(z.score, lower.tail = F)*2
  
  return(pv)
  
  }
  
  
  res.t90$heterog_p.value <- heterog.test.fun(res.t90)
  res.odi$heterog_p.value <- heterog.test.fun(res.odi)
  
  res <- rbind(res.t90,res.odi)

  # Save results 
  fwrite(res, file=paste0(results.folder,"cor_hb_stratified.tsv"))
  
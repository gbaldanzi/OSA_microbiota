# Project sleep apnea and gut microbiota 

# Gabriel Baldanzi 

# Association between T90/ODI and species stratified by the hemoglobin

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


  # Species associated with T90 or ODI in the main model 
  res.main.model <- fread(paste0(results.folder,"cor.bmi_all.var_mgs.tsv"))
  
  mgs.t90 <- res.main.model[q.value<.05 & exposure=="t90",MGS]
  mgs.odi <- res.main.model[q.value<.05 & exposure=="odi",MGS]
  
  
  # Model covariates 
  main.model <- c("age", "Sex", "Alkohol","smokestatus",
                  "plate","shannon","BMI")
  
  
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
  res.t90.hb.low <- lapply(mgs.t90, cor.boot, x = "t90", z=main.model, data=pheno[Hbgroup =="low",])
  res.t90.hb.low <- do.call(rbind, res.t90.hb.low)
  res.t90.hb.low$Hbgroup <- "low"
  
  ## High HB
  res.t90.hb.high <- lapply(mgs.t90, cor.boot, x = "t90", z=main.model, data=pheno[Hbgroup =="high",])
  res.t90.hb.high <- do.call(rbind, res.t90.hb.high)
  res.t90.hb.high$Hbgroup <- "high"
  
  res.t90 <- rbind(res.t90.hb.high, res.t90.hb.low)
  
  res.t90 <- pivot_wider(res.t90, id_cols = c(outcome,exposure), names_from = Hbgroup, 
                           values_from = c(rho,se,conf.int, p.value, covariates))
  
  
  # ODI analysis 
  ## Low HB
  res.odi.hb.low <- lapply(mgs.odi, cor.boot, x = "odi", z=main.model, data=pheno[Hbgroup =="low",])
  res.odi.hb.low <- do.call(rbind, res.odi.hb.low)
  res.odi.hb.low$Hbgroup <- "low"
  
  ## High HB
  res.odi.hb.high <- lapply(mgs.odi, cor.boot, x = "odi", z=main.model, data=pheno[Hbgroup =="high",])
  res.odi.hb.high <- do.call(rbind, res.odi.hb.high)
  res.odi.hb.high$Hbgroup <- "high"
  
  res.odi <- rbind(res.odi.hb.high, res.odi.hb.low)
  
  res.odi <- pivot_wider(res.odi, id_cols = c(outcome,exposure), names_from = Hbgroup, 
                         values_from = c(rho,se,conf.int, p.value, covariates))
  
  # Heterogeneinity test ####
  
  heterog.test.fun <- function(res){
    
  db <- (res$rho_low - res$rho_high)^2
  se <-(res$se_low^2)+(res$se_high^2)
  td <- db/se
  pv <- 1- pchisq(td, df = 1)
  
  return(pv)
  
  }
  
  
  res.t90$heterog_p.value <- heterog.test.fun(res.t90)
  res.odi$heterog_p.value <- heterog.test.fun(res.odi)
  
  # Multiple testing adjustment (FDR)
  res <- rbind(res.t90,res.odi)
  
  res$heterog_q.value <- p.adjust(res$heterog_p.value, method = "BH")
  
  # Save results 
  fwrite(res, file=paste0(results.folder,"cor_hb_stratified.tsv"))
  
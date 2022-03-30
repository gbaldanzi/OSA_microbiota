# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Created: 2022-03-23

# This script will investigate the association between those species that were associated with T90/ODI
# and SBP/DPB/Hb1Ac in SCAPIS-Uppsala participants

results.folder = "/home/baldanzi/Sleep_apnea/Results/"

  # load packages
  library(data.table)
  library(rio)

  # load functions
  source("Script0_functions/Spearman.correlation.function.R")


  # import data ####
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  
  # import the significant species 
  osa.mgs <- readRDS(paste0(results.folder,'mgs.m1.rds'))
  
  # import pathway enrichment analysis results ####
  res.pos <- fread(paste0(results.folder,"ea_GMM_pos.tsv"))
  res.neg <- fread(paste0(results.folder,"ea_GMM_neg.tsv"))
  
  osa.gmm <- unique(res.pos[q.value<.05,pathway], res.neg[q.value<.05,pathway])
  
  
  # Covariates 
  covariates <- c("age", "Sex", "Alkohol","smokestatus",
                  "Fibrer","Energi_kcal" ,"leisurePA", "placebirth",
                  "plate","shannon")
  
  
  # Correlation of significant species with SBP and DBP ####
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = pheno[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.bp <- lapply(outcomes,spearman.function, 
                           x1=c(osa.mgs,osa.gmm),
                           covari = covariates,
                           data = temp.data[hypermed=="no",])
  
  res.bp <- do.call(rbind,res.bp)
  

  # Correlation with HbA1c ####
  # unique(pheno[,Hba1cDetectionLimit])
  # [1] "WITHIN_LIMIT" NA 
  
  outcomes = c("Hba1cFormattedResult")
  
  temp.data = pheno[!is.na(Hba1cFormattedResult),]
  
  res.hb <- lapply(outcomes,spearman.function, 
                   x1=c(osa.mgs,osa.gmm),
                   covari = covariates,
                   data = temp.data[metformin=="no",])
  
  res.hb <- do.call(rbind, res.hb)
  
  res <- rbind(res.bp, res.hb)
  
 
  # Multiple testing adjustment 
  
  setDT(res)
  res[x=="SBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="DBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="Hba1cFormattedResult", q.value := p.adjust(p.value,method = "BH")]
  
  setnames(res, c("x","exposure"), c("MGS_features", "outcomes"))
  
  # Save results ####
  fwrite(res, file=paste0(results.folder, "cor.sig.mgs.gmm_bphb.tsv"))
  

  
  
  # OSA adjusted  ####
  
  ahi.model = c(covariates, "ahi","t90")
  
  # Correlation with SBP and DBP 
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = pheno[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.bp <- lapply(outcomes,spearman.function, 
                   x1=c(osa.mgs,osa.gmm),
                   covari = ahi.model,
                   data = temp.data[hypermed=="no",])
  
  res.bp <- do.call(rbind,res.bp)
  
  # Correlation with HbA1c 
  
  outcomes = c("Hba1cFormattedResult")
  
  temp.data = pheno[!is.na(Hba1cFormattedResult),]
  
  res.hb <- lapply(outcomes,spearman.function, 
                   x1=c(osa.mgs,osa.gmm),
                   covari = ahi.model,
                   data = temp.data[metformin=="no",])
  
  res.hb <- do.call(rbind, res.hb)
  
  res <- rbind(res.bp, res.hb)
  
  # Multiple testing adjustment 
  
  setDT(res)
  res[x=="SBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="DBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="Hba1cFormattedResult", q.value := p.adjust(p.value,method = "BH")]
  
  setnames(res, c("x","exposure"), c("MGS_features", "outcomes"))
  
  # Save results
  fwrite(res, file=paste0(results.folder, "cor.ahi.sig.mgs.gmm_bphb.tsv"))
  
  
  # BMI adjusted ####
  
  covariates.bmi = c(covariates, "BMI", "ahi", "t90")
  
  # Correlation with SBP and DBP
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = pheno[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.bp <- lapply(outcomes,spearman.function, 
                   x1=c(osa.mgs,osa.gmm),
                   covari = covariates.bmi,
                   data = temp.data[hypermed=="no",])
  
  res.bp <- do.call(rbind,res.bp)
  
  # Correlation with HbA1c ####
  
  outcomes = c("Hba1cFormattedResult")
  
  temp.data = pheno[!is.na(Hba1cFormattedResult),]
  
  res.hb <- lapply(outcomes,spearman.function, 
                   x1=c(osa.mgs,osa.gmm),
                   covari = covariates.bmi,
                   data = temp.data[metformin=="no",])
  
  res.hb <- do.call(rbind, res.hb)
  
  res <- rbind(res.bp, res.hb)
  
  
  # Multiple testing adjustment 
  
  setDT(res)
  res[x=="SBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="DBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="Hba1cFormattedResult", q.value := p.adjust(p.value,method = "BH")]
  
  setnames(res, c("x","exposure"), c("MGS_features", "outcomes"))
  
  # Save results ####
  fwrite(res, file=paste0(results.folder, "cor.bmi.sig.mgs.gmm_bphb.tsv"))
  
  
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Created: 2022-03-14

# This script will investigate the association between those species that were associated with T90/ODI
# and SBP/DPB/Hb1Ac. 


# The analysis will be conducted in Uppsala and Malmö participants. 

results.folder = "/home/baldanzi/Sleep_apnea/Results/"

  # load packages
  library(data.table)

  # load functions
  source("Script0_functions/Spearman.correlation.function.R")


  # load data ####
  input <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'
  phenofull <- readRDS(paste0(input,"phenotype_upp_malm.rds"))
  
  # import the significant species 
  osa.mgs <- readRDS(paste0(results.folder,'mgs.m1.rds'))
  
  
  # Covariates 
  covariates <- c("age", "Sex", "Alkohol","smokestatus",
                  "Fibrer","Energi_kcal" ,"leisurePA", "placebirth",
                  "plate","Site" )
  
  
  # Correlation of significant species with SBP and DBP ####
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = phenofull[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.bp <- lapply(outcomes,spearman.function, 
                           x1=osa.mgs,
                           covari = covariates,
                           data = temp.data[hypermed=="no",])
  
  res.bp <- do.call(rbind,res.bp)
  
  # Correlation with HbA1c ####
  # unique(phenofull[,Hba1cDetectionLimit])
  # [1] "WITHIN_LIMIT" NA 
  
  outcomes = c("Hba1cFormattedResult")
  
  temp.data = phenofull[!is.na(Hba1cFormattedResult),]
  
  res.hb <- lapply(outcomes,spearman.function, 
                   x1=osa.mgs,
                   covari = covariates,
                   data = temp.data[metformin=="no",])
  
  res.hb <- do.call(rbind, res.hb)
  
  #Importing enrichment analysis results ####
  res.pos <- fread(paste0(results.folder,"ea_GMM_pos.tsv"))
  res.neg <- fread(paste0(results.folder,"ea_GMM_neg.tsv"))
  
  osa.gmm <- unique(res.pos[q.value<.05,pathway], res.neg[q.value<.05,pathway])
  
  # Correlation of GMMs with SBP and DBP ####
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = phenofull[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.gmm.bp <- lapply(outcomes, spearman.function, 
                   x1 = osa.gmm,
                   covari = covariates,
                   data = temp.data[hypermed=="no",])
  
  res.gmm.bp <- do.call(rbind,res.gmm.bp)
  
  
  # Correlation of GMM with HbA1c ####
  
  outcomes <-  c("Hba1cFormattedResult")
  
  temp.data <-  phenofull[!is.na(Hba1cFormattedResult),]
  
  res.gmm.hb <- lapply(outcomes,spearman.function, 
                   x1 = osa.gmm,
                   covari = covariates,
                   data = temp.data[metformin=="no",])
  
  res.gmm.hb <- do.call(rbind, res.gmm.hb)
  
  res <- rbind(res.bp, res.hb, res.gmm.bp, res.gmm.hb)
  
 
  # Multiple testing adjustment 
  
  setDT(res)
  res[x=="SBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="DBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="Hba1cFormattedResult", q.value := p.adjust(p.value,method = "BH")]
  
  setnames(res, c("x","exposure"), c("MGS_features", "outcomes"))
  
  # Save results ####
  fwrite(res, file=paste0(results.folder, "cor.sig.mgs.gmm_bphb.tsv"))
  

  # Model adjusted for BMI too
  
  covariates.bmi = c(covariates, "BMI")
  
  # Correlation with SBP and DBP ####
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = phenofull[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.bp <- lapply(outcomes,spearman.function, 
                   x1=osa.mgs,
                   covari = covariates.bmi,
                   data = temp.data[hypermed=="no",])
  
  res.bp <- do.call(rbind,res.bp)
  
  # Correlation with HbA1c ####
  
  outcomes = c("Hba1cFormattedResult")
  
  temp.data = phenofull[!is.na(Hba1cFormattedResult),]
  
  res.hb <- lapply(outcomes,spearman.function, 
                   x1=osa.mgs,
                   covari = covariates.bmi,
                   data = temp.data[metformin=="no",])
  
  res.hb <- do.call(rbind, res.hb)
  
  
  # Correlation of GMMs with SBP and DBP ####
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = phenofull[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.gmm.bp <- lapply(outcomes, spearman.function, 
                       x1 = osa.gmm,
                       covari = covariates.bmi,
                       data = temp.data[hypermed=="no",])
  
  res.gmm.bp <- do.call(rbind,res.gmm.bp)
  
  
  # Correlation of GMM with HbA1c ####
  
  outcomes <-  c("Hba1cFormattedResult")
  
  temp.data <-  phenofull[!is.na(Hba1cFormattedResult),]
  
  res.gmm.hb <- lapply(outcomes,spearman.function, 
                       x1 = osa.gmm,
                       covari = covariates.bmi,
                       data = temp.data[metformin=="no",])
  
  res.gmm.hb <- do.call(rbind, res.gmm.hb)
  
  res <- rbind(res.bp, res.hb, res.gmm.bp, res.gmm.hb)
  
  
  # Multiple testing adjustment 
  
  setDT(res)
  res[x=="SBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="DBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res[x=="Hba1cFormattedResult", q.value := p.adjust(p.value,method = "BH")]
  
  setnames(res, c("x","exposure"), c("MGS_features", "outcomes"))
  
  # Save results ####
  fwrite(res, file=paste0(results.folder, "cor.bmi.sig.mgs.gmm_bphb.tsv"))

  
  
  # AHI adjusted  ####
  # Model adjusted for BMI too
  
  ahi.model = c(covariates, "BMI", "ahi","t90")
  ahi.model <- ahi.model[-which(ahi.model=="Site")]
  
  # Correlation with SBP and DBP 
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = phenofull[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.bp <- lapply(outcomes,spearman.function, 
                   x1=c(osa.mgs,osa.gmm),
                   covari = ahi.model,
                   data = temp.data[hypermed=="no",])
  
  res.bp <- do.call(rbind,res.bp)
  
  # Correlation with HbA1c 
  
  outcomes = c("Hba1cFormattedResult")
  
  temp.data = phenofull[!is.na(Hba1cFormattedResult),]
  
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
  
  # Save results ####
  fwrite(res, file=paste0(results.folder, "cor.ahi.sig.mgs.gmm_bphb.tsv"))
  
  
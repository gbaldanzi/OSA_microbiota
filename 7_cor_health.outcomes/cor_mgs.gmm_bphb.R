# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# This script will investigate the association between those species that were 
# associated with T90/ODI and SBP/DPB/Hb1Ac in SCAPIS-Uppsala participants

# It will also investigate the association between the GMM enriched in the T90-species 
# associations and SBP/DPB/Hb1Ac in SCAPIS-Uppsala participants

  library(data.table)

  # Folders 
  results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
  work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'
  

  # load functions
  source("0_functions/Spearman.correlation.function.R")

  # Import data
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))
  

  # Import the significant species 
  osa.mgs <- readRDS(paste0(results.folder,'mgs.m1.rds'))
  
  # Import pathway enrichment analysis results ####
  res.ea.gmm <- fread(paste0(results.folder,"ea_GMM.tsv"))

  osa.gmm <- unique(res.ea.gmm[q.value<.05,pathway])
  
  
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
  
  
  # OSA+BMI adjusted ####
  
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
  
  
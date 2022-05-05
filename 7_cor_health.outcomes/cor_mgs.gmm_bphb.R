# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# This script will investigate the correlation between the species that were 
# associated with T90/ODI and the health outcomes SBP/DPB/Hb1Ac in SCAPIS-Uppsala participants

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
  mgs.fdr <- readRDS(paste0(results.folder,'mgs.m1.rds'))
  
  # Import pathway enrichment analysis results ####
  # res.ea.gmm <- fread(paste0(results.folder,"ea_GMM.tsv"))

  # osa.gmm <- unique(res.ea.gmm[q.value<.05,pathway])
  
  
  # Covariates 
  covariates <- c("age", "Sex", "Alkohol","smokestatus",
                  "Fibrer","Energi_kcal" ,"leisurePA", "placebirth",
                  "plate","t90","ahi","odi")
  
  # In the analysis with SBP and DBP, we excluded participants that self-reported
  # medication use for hypertension 
  # In the analysis with HbA1c, we excluded participants with self-reported 
  # medication use for diabetes 
  
  # Correlation of significant species with SBP and DBP ####
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = pheno[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.bp <- lapply(outcomes,spearman.function, 
                   x1=c(mgs.fdr,osa.gmm),
                   covari = covariates,
                   data = temp.data[hypermed=="no",])
  
  res.bp <- do.call(rbind,res.bp)
  
  # Correlation with HbA1c 
  
  outcomes = c("Hba1cFormattedResult")
  
  temp.data = pheno[!is.na(Hba1cFormattedResult),]
  
  res.hb <- lapply(outcomes,spearman.function, 
                   x1=c(mgs.fdr,osa.gmm),
                   covari = covariates,
                   data = temp.data[diabmed=="no",])
  
  res.hb <- do.call(rbind, res.hb)
  
  res.osa <- rbind(res.bp, res.hb)
  
  # Multiple testing adjustment 
  
  setDT(res.osa)
  res.osa[x=="SBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res.osa[x=="DBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res.osa[x=="Hba1cFormattedResult", q.value := p.adjust(p.value,method = "BH")]
  
  res.osa[,model:="OSA model"]
  
  
  
  # OSA+BMI adjusted ####
  
  covariates.bmi = c(covariates, "BMI")
  
  # Correlation with SBP and DBP
  
  outcomes = c("SBP_Mean", "DBP_Mean")
  
  temp.data = pheno[!is.na(SBP_Mean) & !is.na(DBP_Mean),]
  
  res.bp <- lapply(outcomes,spearman.function, 
                   x1=c(mgs.fdr,osa.gmm),
                   covari = covariates.bmi,
                   data = temp.data[hypermed=="no",])
  
  res.bp <- do.call(rbind,res.bp)
  
  # Correlation with HbA1c ####
  
  outcomes = c("Hba1cFormattedResult")
  
  temp.data = pheno[!is.na(Hba1cFormattedResult),]
  
  res.hb <- lapply(outcomes,spearman.function, 
                   x1=c(mgs.fdr,osa.gmm),
                   covari = covariates.bmi,
                   data = temp.data[diabmed=="no",])
  
  res.hb <- do.call(rbind, res.hb)
  
  res.osa.bmi <- rbind(res.bp, res.hb)
  
  
  # Multiple testing adjustment 
  
  setDT(res.osa.bmi)
  res.osa.bmi[x=="SBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res.osa.bmi[x=="DBP_Mean", q.value := p.adjust(p.value,method = "BH")]
  res.osa.bmi[x=="Hba1cFormattedResult", q.value := p.adjust(p.value,method = "BH")]
  
  res.osa.bmi[,model:="OSA and BMI adjusted"]
  
  
  
  # Save results 
  
  res <- rbind(res.osa, res.osa.bmi)
  
  setnames(res, c("x","exposure"), c("MGS_features", "outcomes"))
  
  fwrite(res, file=paste0(results.folder, "cor_mgs.gmm_bphb.tsv"))
  
  
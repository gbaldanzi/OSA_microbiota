# Project Sleep Apnea and gut microbiota 

# 25 Jan 2022

# Supplementary tables 

  library(xlsx)
  library(tidyverse)
  library(vegan)
  library(data.table)
  library(stringr)
  message("finished loading pckgs")
  rm(list=ls())
  t0 <- Sys.time()
  print(t0)
  if(file.exists('Supp.tables.xlsx')){file.remove('Supp.tables.xlsx')}

  # Function for rounding small values into scientific format 
  round.large <- function(x){
    x[x>=0.001] <- round(x[x>=0.001],3)
    x[x<0.001] <- formatC(x[x<0.001],digits = 2, format = "e")
    return(x)
  }
  
  # Folders
  input <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results'

# Table S3. Association between OSA and alpha-diveristy (Shannon index) ####
  
  
  table.s3 <- fread(paste0(input,"cor_all.var_alpha.tsv"))
  table.s3[,"p-value":=round.large(p.value)]
  table.s3[,"Spearman's correlation":=round(rho,3)]
  table.s3[,exposure := toupper(exposure)]
  
  table.s3 <- table.s3[,c("exposure","Spearman's correlation","p-value","N","model")]
  table.s3 <- table.s3[model %in% c("main.model", "extended.model")]
  table.s3[model =="main.model", model := "Main Model"]
  table.s3[model == "extended.model", model := "Extended Model"]
  
  write.xlsx2(table.s3, "Supp.tables.xlsx", sheetName="Table S3", col.names=T,
              row.names=F, append=F)
  rm(table.s3)

# Table S4. Pairwise comparisons of beta-diversity across Sleep apnea severity categories ####

  list.perma.res <- list.files(path=input, pattern = "pairwise.perma.results")
  
  table.perma.res <- lapply(list.perma.res, function(x){
    table <- readRDS(paste0(input,x))
    table <- table$table
    table$group.1 <-  do.call(rbind,strsplit(table$Comparison,"_", fixed = T))[,1]
    table$group.2 <-  do.call(rbind,strsplit(table$Comparison,"_", fixed = T))[,2]
    table <- table[,c("group.1","group.2","R2 (%)","p-value")]
    return(table)
  })
  
  names(table.perma.res) <- substring(list.perma.res, 24,26)
  
  table.perma.res$ahi$group.1[table.perma.res$ahi$group.1=="noOSA"] <- "No OSA"
  
  table.perma.res$t90$group.1[table.perma.res$t90$group.1=="t0"] <- "T90=0"
  
  white.space <- data.frame(group.1="",group.2="",R2= "", p.value="")
  names(white.space) <- names(table.perma.res$ahi)
  
  table.s4 <- rbind(table.perma.res$ahi, white.space, table.perma.res$t90,
                    white.space, table.perma.res$odi)
  
  write.xlsx2(table.s4, "Supp.tables.xlsx", sheetName="Table S4", col.names=T,
              row.names=F, append=T)
  rm(table.s4)
  rm(table.perma.res)
  rm(list.perma.res)
  rm(white.space)
  

# Table S5 - results from the main model without BMI for all 3 OSA parameters ####

  # Results from main model without BMI 
  res <- fread(paste0(input,"cor_all.var_mgs.tsv"))
  
  #taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  
  #res <- merge(res,taxonomy, by.x="MGS", by.y="maintax_mgs", all.x=T, all.y=F)
  
  res[,MGS:=paste0(MainTax," (",mgs,")")]
  
  res[,exposure:=toupper(exposure)]
  
  res[,rho:=round(rho,3)]
  res[,c("p.value","q.value") := lapply(.SD,round.large) , 
      .SDcols = c("p.value","q.value")]

  
  setnames(res,c("MGS","rho","q.value","p.value"),
           c("Metagenomic species","Spearman's correlation","q-value","p-value"))
  
  var.table <- c("Metagenomic species", "exposure", "Spearman's correlation", "p-value",
                 "q-value", "N","subspecies","species",
                 "genus", "family", "order","class","phylum","superkingdom")
  
  table.res <- res[,var.table,with=F]
  
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Table S5", col.names=T,
              row.names=F, append=T)
  
  # Suppl Table 6 - results from the main model with BMI for all 3 phenotypes  ####
  
  # Results from model2 
  res <- fread(paste0(input,"cor.bmi_all.var_mgs.tsv"))
  
  res[,MGS:=paste0(MainTax," (",mgs,")")]
  
  res[,exposure:=toupper(exposure)]
  
  res[,rho:=round(rho,3)]
  res[,c("p.value","q.value") := lapply(.SD,round.large) , 
      .SDcols = c("p.value","q.value")]
  
  
  setnames(res,c("MGS","rho","q.value","p.value"),
           c("Metagenomic species","Spearman's correlation","q-value","p-value"))
  
  var.table <- c("Metagenomic species", "exposure", "Spearman's correlation", "p-value",
                 "q-value", "N","subspecies","species",
                 "genus", "family", "order","class","phylum","superkingdom")
  
  table.res <- res[,var.table,with=F]
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Table S6", col.names=T,
              row.names=F, append=T)
  
  
  
  # Suppl Table 7 - results from the imputation analysis ####
  message("Results from the imputation analysis")
  print(Sys.time()-t0)
  
  # Results from the AHI-imputed analysis 
  res <- fread(paste0(input,"cor_ahi_imput_mgs.tsv"))
  
  res[,MGS:=gsub("_",".",MGS)]
  
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  
  res <- merge(res,taxonomy, by.x="MGS", by.y="mgs", all.x=T, all.y=F)
  
  res[,MGS:=paste0(MainTax," (",MGS,")")]
  
  res[,exposure:=toupper(exposure)]
  
  res[,rho:=round(rho,3)]
  res[,c("p_value","q_value") := lapply(.SD,round.large) , 
      .SDcols = c("p_value","q_value")]
  
  
  setnames(res,c("MGS","rho","q_value","p_value"),
           c("Metagenomic species","Spearman's correlation","q-value","p-value"))
  
  var.table <- c("Metagenomic species", "exposure", "Spearman's correlation", "p-value",
                 "q-value", "N","subspecies","species",
                 "genus", "family", "order","class","phylum","superkingdom")
  
  table.res <- res[,var.table,with=F]
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Table S7", col.names=T,
              row.names=F, append=T)
  
  
  
  # Suppl Table 8. Sensitivity analyses ####
  message("Table for the Sensitivity anlaysis")
  print(Sys.time()-t0)
  
  # Main model results 
  res.bmi <- fread(paste0(input,"cor.bmi_all.var_mgs.tsv"))
  mgs.t90 <- res.bmi[q.value<.05 & exposure =="t90",MGS]
  mgs.odi <- res.bmi[q.value<.05 & exposure =="odi",MGS]
  
  # Extended model results 
  res <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  
  res1 <- res[MGS %in% mgs.t90 & exposure %in% c("t90")]
  res2 <- res[MGS %in% mgs.odi & exposure %in% c("odi")]
  res <- rbind(res1,res2)
  
  res[,MGS:=paste0(MainTax," (",mgs,")")]
  
  res[,exposure:=toupper(exposure)]
  
  # rounding 
  res[,rho:=round(rho,3)]
  res[,c("p.value","q.value") := lapply(.SD,round.large) , .SDcols = c("p.value","q.value")]
  
  # change names for the final table
  setnames(res,c("rho","p.value","N"), c("extend.model_correlation","extend.model_p-value",
                                         "extend.model_N"))
  
  
  # Medication model results 
  res.med <- fread(paste0(input,"cor.med_all.var_mgs.tsv"))
  
  res1 <- res.med[MGS %in% mgs.t90 & exposure %in% c("t90")]
  res2 <- res.med[MGS %in% mgs.odi & exposure %in% c("odi")]
  res.med <- rbind(res1,res2)
  
  res.med[,MGS:=paste0(MainTax," (",mgs,")")]
  
  res.med[,exposure:=toupper(exposure)]
  
  # rounding 
  res.med[,rho:=round(rho,3)]
  res.med[,c("p.value","q.value") := lapply(.SD,round.large) , .SDcols = c("p.value","q.value")]
  
  # change names for the final table
  res.med <- res.med[,.(MGS,exposure,rho,p.value,N)]
  setnames(res.med,c("rho","p.value","N"), c("medication_correlation","medication_p-value",
                                             "medication_N"))
  
  
  # Antibiotic results 
  res.atb <- fread(paste0(input,"corsaatb_all.var_mgs.tsv"))
  
  res1 <- res.atb[MGS %in% mgs.t90 & exposure %in% c("t90")]
  res2 <- res.atb[MGS %in% mgs.odi & exposure %in% c("odi")]
  res.atb <- rbind(res1,res2)
  
  res.atb <- merge(res.atb,taxonomy, by.x="MGS", by.y="maintax_mgs", all.x=T, all.y=F)
  
  res.atb[,MGS:=paste0(MainTax," (",mgs,")")]
  
  res.atb[,exposure:=toupper(exposure)]
  
  # rounding 
  res.atb[,rho:=round(rho,3)]
  res.atb[,c("p.value","q.value") := lapply(.SD,round.large) , .SDcols = c("p.value","q.value")]
  
  # change names for the final table
  res.atb <- res.atb[,.(MGS,exposure,rho,p.value,N)]
  setnames(res.atb,c("rho","p.value","N"), c("antibiotic_correlation","antibiotic_p-value","antibiotic_N"))
  
  
  # Lung disease results 
  res.lungdisease <- fread(paste0(input,"corsalungdisease_all.var_mgs.tsv"))
  
  res1 <- res.lungdisease[MGS %in% mgs.t90 & exposure %in% c("t90")]
  res2 <- res.lungdisease[MGS %in% mgs.odi & exposure %in% c("odi")]
  res.lungdisease <- rbind(res1,res2)
  
  res.lungdisease <- merge(res.lungdisease,taxonomy, by.x="MGS", by.y="maintax_mgs", all.x=T, all.y=F)
  
  res.lungdisease[,MGS:=paste0(MainTax," (",mgs,")")]
  
  res.lungdisease[,exposure:=toupper(exposure)]
  
  # rounding 
  res.lungdisease[,rho:=round(rho,3)]
  res.lungdisease[,c("p.value","q.value") := lapply(.SD,round.large) , .SDcols = c("p.value","q.value")]
  
  # change names for the final table
  res.lungdisease <- res.lungdisease[,.(MGS,exposure,rho,p.value,N)]
  setnames(res.lungdisease,c("rho","p.value","N"), c("lung_correlation","lung_p-value","lung_N"))
  
  
  
  # Merge all sensitivity analysis results 
  res <- merge(res,res.med, by =c("MGS","exposure"))
  res <- merge(res,res.atb, by = c("MGS","exposure"))
  
  
  # selecting variables for the final table
  setnames(res,"MGS","Metagenomic species")
  var.table <- c("Metagenomic species", "exposure", "extend.model_correlation", 
                 "extend.model_p-value","extend.model_N","medication_correlation",
                 "medication_p-value","medication_N", "antibiotic_correlation","antibiotic_p-value",
                 "antibiotic_N","lung_correlation","lung_p-value",
                 "lung_N",
                 "species","genus", "family", "order","class","phylum","superkingdom")
  
  table.res <- res[,var.table,with=F]
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Table S8", col.names=T,
              row.names=F, append=T)
  rm(table.res)
  rm(res.atb)
  rm(res.med)
  rm(res.lungdisease)
  rm(res)
  
  
  
  # Suppl Table 9 - GMM enrichment analysis ####
  message("Table for the GMM enrichment analysis results")
  
  
  res <- fread(paste0(input,"ea_GMM.tsv"))
  
  res[,exposure:=toupper(exposure)]
  
  
  res <- res[,.(exposure,correlation,pathway,Name,HL1,HL2,pval,q.value,NES,size)]
  
  res[,c("pval","q.value","NES"):=lapply(.SD,round.large), .SDcols=c("pval","q.value","NES")]
  
  setnames(res,c("pathway","pval","q.value"),c("Gut metabolic module","p-value","q-value"))
  
  setDF(res)
  write.xlsx2(res, "Supp.tables_GMM.xlsx", sheetName="Table S9", col.names=T,
              row.names=F, append=F)

  
  message(paste("File name = Supp.tables.xlsx, saved at",getwd()))
  message(paste("File name = Supp.tables_GMM.xlsx, saved at",getwd()))
  
  
  # Suppl Table 10 - GMM-metabolites Spearman correlation ####
  
  cor_gmm_metabolites <- fread(paste0(input,"cor_gmm_metabolites.tsv"))
  
  gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  gmm.names[,Name:=str_to_title(Name)]
  gmm.names[,Name:=gsub("Ii","II",Name)]
  
  cor_gmm_metabolites <- merge(cor_gmm_metabolites, gmm.names[,.(Module,Name)],
                               by.x= "modules", by.y = "Module", all.x=T, all.y=F)
  
  input2 <- "/proj/sens2019512/SCAPIS_org/SCAPIS/metabolon_final/clean/"
  annotation <- fread(paste0(input2, 'scapis_merged_annotations_batchnorm_clean.tsv'))
  
  cor_gmm_metabolites <- merge(cor_gmm_metabolites, annotation[,.(MET_ID, CHEMICAL_NAME,SUB_PATHWAY)], 
               by.x = "metabolites", by.y = "MET_ID", all.x=T, all.y=F)
  
  res <- cor_gmm_metabolites[,.(modules, Name,CHEMICAL_NAME,SUB_PATHWAY,rho,p.value,q.value,N)]
  
  res[!is.na(p.value), p.value := round.large(p.value)]
  res[!is.na(q.value), q.value := round.large(q.value)]
  
  setnames(res, c("modules", "CHEMICAL_NAME","SUB_PATHWAY", "rho", "p.value", "q.value"),
           c("Gut metabolic module", "Metabolite","Metabolite group","Spearman's correlation",
             "p-value","q-value" ))
  
  write.table(res, "Supp.tables_S10.tsv", sep = "\t", row.names = F)
  
  message(paste("File name = Supp.tables_S10.tsv, saved at",getwd()))
  
  # Suppl Table 11 - GMM-metabolites enrichment ####
  
  ea_gmm_subclass <- fread(paste0(input,"ea_gmm_subclass.tsv"))
  
  ea_gmm_subclass <- merge(ea_gmm_subclass, gmm.names[,.(Module,Name)],
                               by.x= "modules", by.y = "Module", all.x=T, all.y=F)
  
  fix.subclass.name.fun <- function(char){
    char <- gsub("__PC_", " (PC)",char)
    char <- gsub("___","_",char)
    char <- gsub("__","_",char)
    char <- gsub("Drug", "Drug -", char)
    char <- gsub("_"," ",char)
    char <- gsub(" PI ", " (PI)",char)
    char <- gsub(" PE ", " (PE)",char)
    char <- gsub("Analgesics","Analgesics,", char)
    return(char)
  }
  
  res <- ea_gmm_subclass[,.(modules,Name,direction,subclass,estimate,p.value,q.value,size)]
  
  res[, subclass := fix.subclass.name.fun(subclass)]
  
  res[,c("p.value","q.value","estimate"):=lapply(.SD,round.large), .SDcols=c("p.value","q.value","estimate")]
  
  setnames(res,c("modules","subclass","estimate","p.value","q.value"),
           c("Gut metabolic module","Metabolites group","NES","p-value","q-value"))

  write.xlsx2(res, "Supp.tables_part2.xlsx", sheetName="Table S11", col.names=T,
              row.names=F, append=F)
  
  # Suppl Table 12 - health outcomes ####
  
  res.mgs.bp <- fread(paste0(input, 'cor.sig.mgs.gmm_bphb.tsv'))
  res.mgs.bp[,model:="basic model"]
  res.mgs.bp.ahi <- fread(paste0(input, 'cor.ahi.sig.mgs.gmm_bphb.tsv'))
  res.mgs.bp.ahi[,model:="OSA adjusted"]
  res.mgs.bp.bmi <- fread(paste0(input, 'cor.bmi.sig.mgs.gmm_bphb.tsv'))
  res.mgs.bp.bmi[,model:="OSA and BMI adjusted"]
  
  res <- rbind(res.mgs.bp, res.mgs.bp.ahi, res.mgs.bp.bmi)
  
  res1 <- merge(res[!grep("MF",MGS_features),], taxonomy, by.x="MGS_features", by.y="maintax_mgs",
               all.x=T, all.y=F)
  res2 <- merge(res[grep("MF",MGS_features),], gmm.names[,.(Module,Name)], by.x = "MGS_features", by.y="Module",
               all.x=T, all.y=F)
  
  res1[, microbiota := paste0(MainTax," (",mgs,")")]
  res2[, microbiota := Name]
  
  res1 <- res1[,.(microbiota, outcomes, rho, p.value, q.value, N, model)]
  res2 <- res2[,.(microbiota, outcomes, rho, p.value, q.value, N, model)]
  
  res <- rbind(res1,res2)
  
  res[,outcomes:= factor(outcomes, levels = c("SBP_Mean","DBP_Mean","Hba1cFormattedResult"),
                         labels = c("SBP", "DBP", "HbA1c"))]
  
  res[,c("Spearman's correlation","p-value","q-value") := lapply(.SD, round.large), 
      .SDcols = c("rho","p.value","q.value")]
  
  res <- res[,c("microbiota","outcomes","Spearman's correlation","p-value",
                "q-value","N","model")]
    
  
  
  write.xlsx2(res, "Supp.tables_part2.xlsx", sheetName="Table S12", col.names=T,
              row.names=F, append=T)
  
  message(paste("File name = Supp.tables_part2.xlsx, saved at",getwd()))
  
  print(Sys.time())
  
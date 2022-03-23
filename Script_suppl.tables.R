# Project Sleep Apnea and gut microbiota 

# 25 Jan 2022

# Supplementary tables 

library(xlsx)
library(tidyverse)
library(vegan)
library(data.table)
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

# Table S3. Association between OSA and alpha-diveristy (Shannon index) ####
  input <- "/home/baldanzi/Sleep_apnea/Results/"
  
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

# Table S4. Pairwise comparisons of beta-diversity across Sleep apnea severity categories

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
  
  
  # Medical model results 
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
  
  
  
  # Merge all sensitivity analysis results 
  res <- merge(res,res.med, by =c("MGS","exposure"))
  res <- merge(res,res.atb, by = c("MGS","exposure"))
  
  
  # selecting variables for the final table
  setnames(res,"MGS","Metagenomic species")
  var.table <- c("Metagenomic species", "exposure", "extend.model_correlation", 
                 "extend.model_p-value","extend.model_N","medication_correlation",
                 "medication_p-value","medication_N", "antibiotic_correlation","antibiotic_p-value",
                 "antibiotic_N","species","genus", "family", "order","class","phylum","superkingdom")
  
  table.res <- res[,var.table,with=F]
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Table S8", col.names=T,
              row.names=F, append=T)
  rm(table.res)
  rm(res.atb)
  rm(res.med)
  rm(res)
  
  
  
  # Suppl Table 9 - GMM enrichment analysis ####
  message("Table for the GMM enrichment analysis results")
  
  
  res.pos <- fread(paste0(input,"ea_GMM_pos.tsv"))
  res.neg <- fread(paste0(input,"ea_GMM_neg.tsv"))
  res.pos[,correlation:="positive"]
  res.neg[,correlation:="negative"]
  
  res <- rbind(res.pos,res.neg)
  
  res[,exposure:=toupper(exposure)]
  
  
  res <- res[,.(exposure,correlation,pathway,Name,HL1,HL2,pval,q.value,NES,size)]
  
  res[,c("pval","q.value","NES"):=lapply(.SD,round.large), .SDcols=c("pval","q.value","NES")]
  
  setnames(res,c("pathway","pval","q.value"),c("Gut metabolic module","p-value","q-value"))
  
  setDF(res)
  write.xlsx2(res, "Supp.tables_GMM.xlsx", sheetName="Table S9", col.names=T,
              row.names=F, append=F)

  
  message(paste("File name = Supp.tables.xlsx, saved at",getwd()))
  message(paste("File name = Supp.tables_GMM.xlsx, saved at",getwd()))
  print(Sys.time())
  
# Project Sleep Apnea 

# Supplementary tables 

library(xlsx)
library(tidyverse)
library(vegan)
library(data.table)

# Table S1. Pairwise comparisons of beta-diversity across Sleep apnea severity categories
  input <- "/home/baldanzi/Sleep_apnea/Results/"
  
  list.res.ahi <- readRDS(paste0(input,"pairwise.perma.results.rds"))
  list.res.ahi <- list.res.ahi[1:6]
  
  list.res.t90 <- readRDS(file=paste0(input,"pairwise.perma.results_t90.rds"))
  
  pairwise.table.fun <- function(lis) {
    
    table <- lapply(lis, function(x){
      x[["variables"]] <-  gsub("no_OSA","noOSA",x[["variables"]])
      l <- x[1,c("variables","Pr(>F)")] %>% 
        separate(variables, c("group.1", "group.2")) %>%
        rename(p.value = "Pr(>F)")
      return(l)
    })
    
    table <- do.call(rbind,table)
    
    return(table)
  }
  
  table.res.ahi <- pairwise.table.fun(list.res.ahi)
  table.res.ahi$group.1[table.res.ahi$group.1=="noOSA"] <- "No OSA"
  
  table.res.t90 <- pairwise.table.fun(list.res.t90)
  table.res.t90$group.1[table.res.t90$group.1=="t0"] <- "T90=0"
  
  
  white.space <- data.frame(group.1="",group.2="",p.value="")
  
  table.s1 <- rbind(table.res.ahi, white.space, table.res.t90)
  
  write.xlsx2(table.s1, "Supp.tables.xlsx", sheetName="Table S1", col.names=T,
              row.names=F, append=F)
  

# Suppl Table 2 - results from the basic model for all 4 phenotypes ####

  # Results from model1 
  res <- fread(paste0(input,"cor_all.var_mgs.tsv"))
  
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  
  res <- merge(res,taxonomy, by.x="MGS", by.y="maintax_mgs", all.x=T, all.y=F)
  
  res[,MGS:=paste0(MainTax," (",mgs,")")]
  
  res[,exposure:=toupper(exposure)]
  
  res[,c("cor.coefficient","p.value","q.value") := lapply(.SD,round,3) , 
      .SDcols = c("cor.coefficient","p.value","q.value")]

  
  setnames(res,c("MGS","cor.coefficient","q.value","p.value"),
           c("Metagenomic species","Spearman's correlation","q-value","p-value"))
  
  var.table <- c("Metagenomic species", "exposure", "Spearman's correlation", "p-value",
                 "q-value", "N","subspecies","species",
                 "genus", "family", "order","class","phylum","superkingdom")
  
  table.res <- res[,var.table,with=F]
  
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Table S2", col.names=T,
              row.names=F, append=T)
  
  # Suppl Table 3 - results from the fully adjusted model for all 3 phenotypes  ####
  
  # Results from model2 
  res <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  
  res[,MGS:=paste0(MainTax," (",mgs,")")]
  
  res[,exposure:=toupper(exposure)]
  
  res[,c("cor.coefficient","p.value","q.value") := lapply(.SD,round,3) , 
      .SDcols = c("cor.coefficient","p.value","q.value")]
  
  
  setnames(res,c("MGS","cor.coefficient","q.value","p.value"),
           c("Metagenomic species","Spearman's correlation","q-value","p-value"))
  
  var.table <- c("Metagenomic species", "exposure", "Spearman's correlation", "p-value",
                 "q-value", "N","subspecies","species",
                 "genus", "family", "order","class","phylum","superkingdom")
  
  table.res <- res[,var.table,with=F]
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Table S3", col.names=T,
              row.names=F, append=T)
  
  # Suppl Table 4. List of Signatures species / Metagenomics information ####
  
  mgs.m2 <- readRDS(paste0(input,'mgs.m2.rds'))
  
  mgs.m2 <- unique(do.call(c,mgs.m2))
  
  
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  prev <- pheno[,mgs.m2,with=F]
  prev[,(mgs.m2) := lapply(.SD,decostand,method="pa"), .SDcols=mgs.m2]
  prev <- data.frame(maintax_mgs=mgs.m2,
                     prevalence=round(t(prev[, lapply(.SD,function(x) sum(x)/nrow(prev)), 
                                             .SDcols=mgs.m2]),3)*100)
  
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  
  taxonomy <- taxonomy[maintax_mgs %in% mgs.m2,]
  
  taxonomy <- merge(taxonomy,prev,by="maintax_mgs")
  
  taxonomy[,MainTax:= paste0(MainTax, " (", mgs,")")]
  
  a = c("mgs","maintax_mgs","Level")
  
  taxonomy <- taxonomy[,-a,with=F]
  
  setnames(taxonomy,"MainTax","OSA signature species")
  
  setcolorder(taxonomy,c("OSA signature species","prevalence"))
  
  
  write.xlsx2(taxonomy, "Supp.tables.xlsx", sheetName="Table S4", col.names=T,
              row.names=F, append=T)
  
  # Suppl. Table 5 - medication use ####
  
  message("Medication use")
  
  res_full <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  setnames(res_full,"cor.coefficient","rho")
  res_sa_med <- fread(paste0(input,"corsa_all.var_mgs.tsv"))
  setnames(res_sa_med,"cor.coefficient","rho")
  mgs.m2 <- readRDS(paste0(input,'mgs.m2.rds'))
  
  names(mgs.m2) <- c("ahi","t90","odi")
  
  round.large <- function(x){
    x[x>=0.001] <- round(x[x>=0.001],3)
    x[x<0.001] <- formatC(x[x<0.001],digits = 2, format = "e")
    return(x)
  }
  
  prepare.sa.table.fun <- function(expo,res,res_sa) { 
  
    stopifnot("'expo' needs to be of class 'character'"=class(expo) %in% "character")  
    
    res_sa1 <- res_sa[exposure==expo & MGS %in% mgs.m2[[expo]],.(MGS,exposure,rho,p.value,N)]
    res1 <- res[exposure==expo & MGS %in% mgs.m2[[expo]],.(MGS,rho,p.value,N,MainTax,mgs)]
    setnames(res1,c('rho','p.value','N'),paste0(c('rho','p.value','N'),"_full_model"))
    res_sa1 <- merge(res_sa1,res1,by="MGS",all.x=T)
    res_sa1[ , MGS:=paste0(MainTax, " (", mgs,")")]
    res_sa1 <- res_sa1[,.(MGS,exposure,rho_full_model,p.value_full_model,N_full_model,rho,p.value,N)]
    res_sa1[,exposure:=toupper(exposure)]
      cols <- c("rho_full_model","rho")
    res_sa1[,(cols):=lapply(.SD,round,digits=3),.SDcols=cols]
      cols <- c("p.value_full_model","p.value")
    res_sa1[,(cols):=lapply(.SD,round.large),.SDcols=cols]
    
    return(res_sa1)
  
  }
  
  
  t <- lapply(c("ahi","t90","odi"), prepare.sa.table.fun
                , res = res_full, res_sa=res_sa_med)
  
  break.fun <- function(expo) {
  
        data.table(MGS=c("",paste0("Signatures species associated to ",expo)),
                       exposure=c("","exposure"),
                       rho_full_model = c("Fully adjusted model","Spearman's correlation"),
                       p.value_full_model = c("Fully adjusted model","p-value"),
                       N_full_model = c("Fullly adjusted model","N"),
                       rho = c("Removed medication users","Spearman's correlation"),
                       p.value = c("Removed medication users","p-value"),
                       N= c("Removed medication users", "N"))
  }
  
  brks <- lapply(c("AHI","T90","ODI"), break.fun)
  
  res <- rbind(brks[[1]], t[[1]], brks[[2]], t[[2]], brks[[3]], t[[3]])
  
  
  write.xlsx2(res, "Supp.tables.xlsx", sheetName="Table S5", col.names=F,
              row.names=F, append=T)
  
  
  
  # Suppl. Table 6 - antibiotic use ####
  
  message("Antibiotic use")
  
  res_full <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  setnames(res_full,"cor.coefficient","rho")
  res_sa_atb <- fread(paste0(input,"cor_sa_atb6m_step2_all.var_mgs.tsv"))
  
  t <- lapply(c("ahi","t90","odi"), prepare.sa.table.fun
              , res = res_full, res_sa=res_sa_atb)
  
  
  break.fun <- function(expo) {
    
    data.table(MGS=c("",paste0("Signatures species associated to ",expo)),
               exposure=c("","exposure"),
               rho_full_model = c("Fully adjusted model","Spearman's correlation"),
               p.value_full_model = c("Fully adjusted model","p-value"),
               N_full_model = c("Fullly adjusted model","N"),
               rho = c("No antibiotic use last 6 months","Spearman's correlation"),
               p.value = c("No antibiotic use last 6 months","p-value"),
               N= c("No antibiotic use last 6 months", "N"))
  }
  
  brks <- lapply(c("AHI","T90","ODI"), break.fun)
  
  res <- rbind(brks[[1]], t[[1]], brks[[2]], t[[2]], brks[[3]], t[[3]])
  
  
  write.xlsx2(res, "Supp.tables.xlsx", sheetName="Table S6", col.names=F,
              row.names=F, append=T)
  
  
  
  # Suppl Table 7 - GMM enrichment analysis ####
  
  res.pos <- fread(paste0(input,"ea_GMM_pos.tsv"))
  res.neg <- fread(paste0(input,"ea_GMM_neg.tsv"))
  res.pos[,correlation:="positive"]
  res.neg[,correlation:="negative"]
  
  res <- rbind(res.pos,res.neg)
  
  gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  res <- merge(res,gmm.names,by.x="pathway", by.y="Module", all.x=T)
  res[,exposure:=toupper(exposure)]
  
  
  res <- res[,.(exposure,correlation,pathway,Name,HL1,HL2,pval,q.value,NES,size)]
  
  res[,c("pval","q.value","NES"):=lapply(.SD,round.large), .SDcols=c("pval","q.value","NES")]
  
  setnames(res,c("pathway","pval","q.value"),c("Gut metabolic module","p-value","q-value"))
  
  write.xlsx2(res, "Supp.tables.xlsx", sheetName="Table S7", col.names=T,
              row.names=F, append=T)
  
  # Suppl Table 8 - GMM enrichment analysis - sensitivity analysis ####
  
  message("GMM enrichment - SA")
  
  res.pos <- fread(paste0(input,"ea_GMM_pos.tsv"))
  cols <- c("NES","pval","q.value")
  res.pos[,(cols):=lapply(.SD,round.large),.SDcols=cols]
  
  bmi.pw <- res.pos[exposure=="BMI" & q.value<.05,pathway]
  
  res.pos_sa <- fread("/home/baldanzi/Sleep_apnea/Results/ea_GMM_pos_sa_metformin.tsv")
  res.pos_sa[,exposure:=toupper(exposure)]
  res.pos_sa[,(cols):=lapply(.SD,round.large),.SDcols=cols]
  

  t90.pw <- res.pos[exposure=="t90" & q.value<.05 & !pathway %in% bmi.pw,
                    .(pathway,pval,q.value,NES)]
  
  names(t90.pw) <- paste0(names(t90.pw), "_main")
  
  t90.pw_sa <- res.pos_sa[exposure=="T90" & pathway %in% t90.pw$pathway,
                          .(pathway,exposure,pval,q.value,NES)]
  
  t90.pw_sa <- merge(t90.pw_sa, t90.pw, by.x="pathway",by.y="pathway_main", all.x=T)
  
  t90.pw_sa <- merge(t90.pw_sa, gmm.names, by.x="pathway",by.y="Module", all.x=T)
  
  #setnames(t90.pw_sa, c("pathway","pval","q.value"), 
           #c("Gut metabolic module","p-value","q-value"))
  
  setcolorder(t90.pw_sa,c("exposure","pathway","Name", "HL1", "HL2",
                          "pval_main","q.value_main","NES_main"))
  
  
  top <- data.table(exposure=c("","exposure"),
                       pathway=c("","Gut metabolic module"),
                       Name=c("","Name"),
                       HL1=c("","HL1"),
                       HL2=c("","HL2"),
                       pval_main = c("Main analysis","p-value"),
                      q.value_main = c("Main analysis","q-value"),
                       NES_main = c("Main analysis","NES"),
                       pval = c("Removed metformin users","p-value"),
                       q.value = c("Removed metformin users","q-value"),
                       NES = c("Removed metformin users","NES"))
  
  

  res <- rbind(top,t90.pw_sa)
  
  message("Saving final version")
  
  
  write.xlsx2(res, "Supp.tables.xlsx", sheetName="Table S8", col.names=F,
              row.names=F, append=T)
  
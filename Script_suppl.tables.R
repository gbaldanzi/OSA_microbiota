# Project Sleep Apnea 

# Supplementary tables 

library(xlsx)

# Table S1. Pairwise comparisons of beta-diversity across Sleep apnea severity categories
  input <- "/home/baldanzi/Sleep_apnea/Results/"
  list.res <- readRDS(paste0(input,"pairwise.perma.results.rds"))
  
  Severe <- c(list.res$no_OSA_Severe[1,"Pr(>F)"],
              list.res$Mild_Severe[1,"Pr(>F)"],
              list.res$Moderate_Severe[1,"Pr(>F)"])
  
  Moderate <- c(list.res$no_OSA_Moderate[1,"Pr(>F)"],
                list.res$Mild_Moderate[1,"Pr(>F)"],"")
  
  Mild <- c(list.res$no_OSA_Mild[1,"Pr(>F)"],"","")
              
  
  table.res <- data.frame(categories=c("No Sleep Apnea", "Mild", "Moderate"),
                          Severe = Severe,
                          Moderate = Moderate,
                          Mild = Mild )
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Pairwise comparisons BC", col.names=T,
              row.names=F, append=F)
  
# Suppl Table 2 - results from the basic model for all 3 phenotypes 

  # Results from model1 
  res <- fread(paste0(input,"cor_all.var_mgs.tsv"))
  
  res[,MGS:=gsub("____"," (",MGS)]
  res[,MGS:=paste0(MGS,")")]
  res[,MGS:=gsub("_"," ",MGS)]
  
  res[,exposure:=toupper(exposure)]
  
  setnames(res,c("MGS","cor.coefficient","q.value"),c("Metagenomic species","rho","FDR_p"))
  
  var.table <- c("Metagenomic species", "exposure", "rho", "p.value", "N","subspecies","species",
                 "genus", "family", "order","class","phylum","superkingdom")
  
  table.res <- res[,var.table,with=F]
  
  setDT(table.res)
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Results Basic Model", col.names=T,
              row.names=F, append=T)
  
  # Suppl Table 3 - results from the fully adjusted model for all 3 phenotypes 
  
  # Results from model2 
  res <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  
  res[,MGS:=gsub("____"," (",MGS)]
  res[,MGS:=paste0(MGS,")")]
  res[,MGS:=gsub("_"," ",MGS)]
  
  res[,exposure:=toupper(exposure)]
  
  setnames(res,c("MGS","cor.coefficient","q.value"),c("Metagenomic species","rho","FDR_p"))
  
  var.table <- c("Metagenomic species", "exposure", "rho", "p.value", "N","subspecies","species",
                 "genus", "family", "order","class","phylum","superkingdom")
  
  table.res <- res[,var.table,with=F]
  
  setDT(table.res)
  
  write.xlsx2(table.res, "Supp.tables.xlsx", sheetName="Results Fully Adjusted Model", col.names=T,
              row.names=F, append=T)
  
  # Suppl Table 4. List of Signatures species / Metagenomics information
  
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  
  mgs.m2 <- readRDS(paste0(input,'mgs.m2.rds'))
  
  mgs.m2 <- unique(c(mgs.m2[[1]],mgs.m2[[2]]))
  
  taxonomy <- taxonomy[maintax_mgs %in% mgs.m2,]
  
  taxonomy[,maintax_mgs:=gsub("____"," (",maintax_mgs)]
  taxonomy[,maintax_mgs:=paste0(maintax_mgs,")")]
  taxonomy[,maintax_mgs:=gsub("_"," ",maintax_mgs)]
  
  a = c("mgs","MainTax","Level")
  
  taxonomy <- taxonomy[,-a,with=F]
  
  setnames(taxonomy,"maintax_mgs","Metagenomics species")
  
  
  write.xlsx2(taxonomy, "Supp.tables.xlsx", sheetName="Taxonomy Signature Species", col.names=T,
              row.names=F, append=T)
  
  # Suppl. Table 5 - medication use ####
  
  message("Medication use")
  
  res <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  setnames(res,"cor.coefficient","rho")
  res_sa <- fread(paste0(input,"corsa_all.var_mgs.tsv"))
  setnames(res_sa,"cor.coefficient","rho")
  mgs.m2 <- readRDS(paste0(input,'mgs.m2.rds'))
  
  round.large <- function(x){
    x[x>=0.001] <- round(x[x>=0.001],3)
    x[x<0.001] <- formatC(x[x<0.001],digits = 2, format = "e")
    return(x)
  }
  
  res_sa1 <- res_sa[exposure=="ahi" & MGS %in% mgs.m2$mgs.fdr.ahi,.(MGS,exposure,rho,p.value,N)]
  res1 <- res[exposure=="ahi" & MGS %in% mgs.m2$mgs.fdr.ahi,.(MGS,rho,p.value)]
  setnames(res1,names(res1[,2:3]),paste0(names(res1[,2:3]),"_full_model"))
  res_sa1 <- merge(res_sa1,res1,by="MGS",all.x=T)
  res_sa1 <- res_sa1[,.(MGS,exposure,rho_full_model,p.value_full_model,rho,p.value,N)]
  res_sa1[,exposure:=toupper(exposure)]
    cols <- c("rho_full_model","rho")
  res_sa1[,(cols):=lapply(.SD,round,digits=3),.SDcols=cols]
    cols <- c("p.value_full_model","p.value")
  res_sa1[,(cols):=lapply(.SD,round.large),.SDcols=cols]
  
  res_sa2 <- res_sa[exposure=="t90" & MGS %in% mgs.m2$mgs.fdr.t90,.(MGS,exposure,rho,p.value,N)]
  res2 <- res[exposure=="t90" & MGS %in% mgs.m2$mgs.fdr.t90,.(MGS,rho,p.value)]
  setnames(res2,names(res2[,2:3]),paste0(names(res2[,2:3]),"_full_model"))
  res_sa2 <- merge(res_sa2,res2,by="MGS",all.x=T)
  res_sa2 <- res_sa2[,.(MGS,exposure,rho_full_model,p.value_full_model,rho,p.value,N)]
  res_sa2[,exposure:=toupper(exposure)]
  cols <- c("rho_full_model","rho")
  res_sa2[,(cols):=lapply(.SD,round,digits=3),.SDcols=cols]
  cols <- c("p.value_full_model","p.value")
  res_sa2[,(cols):=lapply(.SD,round.large),.SDcols=cols]
  
  break1 <- data.table(MGS="Signatures species associated to AHI",
                       exposure="",
                       rho_full_model = "Fully adjusted model",
                       p.value_full_model = "Fully adjustd model",
                       rho = "Removed medication users",
                       p.value = "Removed medication users",
                       N= "Removed medication users")
  
  break2 <- data.table(MGS="Signatures species associated to T90",
                       exposure="",
                       rho_full_model = "Fully adjusted model",
                       p.value_full_model = "Fully adjustd model",
                       rho = "Removed medication users",
                       p.value = "Removed medication users",
                       N= "Removed medication users")
  
  res <- rbind(break1, res_sa1, break2, res_sa2)
  
  res[,MGS:=gsub("____"," (",MGS)]
  res[,MGS:=paste0(MGS,")")]
  res[,MGS:=gsub("_"," ",MGS)]
  setnames(res,"MGS","Metagenomic species")
  
  setDT(res)
  
  write.xlsx2(res, "Supp.tables.xlsx", sheetName="Sens. Analysis - Medications", col.names=T,
              row.names=F, append=T)
  
  
  
  # Suppl. Table 6 - antibiotic use ####
  
  message("Antibiotic use")
  
  res <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  setnames(res,"cor.coefficient","rho")
  res_sa <- fread(paste0(input,"cor_sa_atb6m_step2_all.var_mgs.tsv"))
  
  mgs.m2 <- readRDS(paste0(input,'mgs.m2.rds'))
  

  res_sa1 <- res_sa[exposure=="ahi" & MGS %in% mgs.m2$mgs.fdr.ahi,.(MGS,exposure,rho,p.value,N)]
  res1 <- res[exposure=="ahi" & MGS %in% mgs.m2$mgs.fdr.ahi,.(MGS,rho,p.value)]
  setnames(res1,names(res1[,2:3]),paste0(names(res1[,2:3]),"_full_model"))
  res_sa1 <- merge(res_sa1,res1,by="MGS",all.x=T)
  res_sa1 <- res_sa1[,.(MGS,exposure,rho_full_model,p.value_full_model,rho,p.value,N)]
  res_sa1[,exposure:=toupper(exposure)]
  cols <- c("rho_full_model","rho")
  res_sa1[,(cols):=lapply(.SD,round,digits=3),.SDcols=cols]
  cols <- c("p.value_full_model","p.value")
  res_sa1[,(cols):=lapply(.SD,round.large),.SDcols=cols]
  
  res_sa2 <- res_sa[exposure=="t90" & MGS %in% mgs.m2$mgs.fdr.t90,.(MGS,exposure,rho,p.value,N)]
  res2 <- res[exposure=="t90" & MGS %in% mgs.m2$mgs.fdr.t90,.(MGS,rho,p.value)]
  setnames(res2,names(res2[,2:3]),paste0(names(res2[,2:3]),"_full_model"))
  res_sa2 <- merge(res_sa2,res2,by="MGS",all.x=T)
  res_sa2 <- res_sa2[,.(MGS,exposure,rho_full_model,p.value_full_model,rho,p.value,N)]
  res_sa2[,exposure:=toupper(exposure)]
  cols <- c("rho_full_model","rho")
  res_sa2[,(cols):=lapply(.SD,round,digits=3),.SDcols=cols]
  cols <- c("p.value_full_model","p.value")
  res_sa2[,(cols):=lapply(.SD,round.large),.SDcols=cols]
  
  break1 <- data.table(MGS="Signatures species associated to AHI",
                       exposure="",
                       rho_full_model = "Fully adjusted model",
                       p.value_full_model = "Fully adjustd model",
                       rho = "No antibiotic use last 6 months",
                       p.value = "No antibiotic use last 6 months",
                       N= "No antibiotic use last 6 months")
  
  break2 <- data.table(MGS="Signatures species associated to T90",
                       exposure="",
                       rho_full_model = "Fully adjusted model",
                       p.value_full_model = "Fully adjustd model",
                       rho = "No antibiotic use last 6 months",
                       p.value = "No antibiotic use last 6 months",
                       N= "No antibiotic use last 6 months")
  
  res <- rbind(break1, res_sa1, break2, res_sa2)
  
  setDT(res)
  
  res[,MGS:=gsub("____"," (",MGS)]
  res[,MGS:=paste0(MGS,")")]
  res[,MGS:=gsub("_"," ",MGS)]
  setnames(res,"MGS","Metagenomic species")
  
  write.xlsx2(res, "Supp.tables.xlsx", sheetName="Sens. Analysis - Antibiotic", col.names=T,
              row.names=F, append=T)
  
  # Suppl. Table 7 - Sub-pathways enrichment analysis ####
  
  res.pos <- fread(paste0(input,"ea_subpathways_pos.tsv"))
  res.pos[,correlation:="positive"]
  res.neg <- fread(paste0(input,"ea_subpathways_neg.tsv"))
  res.neg[,correlation:="negative"]
  
  res <- rbind(res.pos,res.neg)
  
  
  res <- res[,.(MGS,correlation,pathway,pval,padj,ES,NES,size)]
  
  res[,MGS:=gsub("____"," (",MGS)]
  res[,MGS:=paste0(MGS,")")]
  res[,MGS:=gsub("_"," ",MGS)]
  res[,MGS:=gsub("42 11","42-11",MGS)]
  
  setnames(res,c("MGS","pathway","pval","padj"),c("Metagenomic species","sub-pathways",
                                                  "p.value","FDR_p"))
  
  setDF(res)
  
  write.xlsx2(res, "Supp.tables.xlsx", sheetName="Sub-pathway enrichment analysis", col.names=T,
              row.names=F, append=T)
  
  
  # Suppl Table 8 - GMM enrichment analysis ####
  
  res.pos <- fread(paste0(input,"ea_GMM_pos.tsv"))
  res.neg <- fread(paste0(input,"ea_GMM_neg.tsv"))
  res.pos[,correlation:="positive"]
  res.neg[,correlation:="negative"]
  
  res <- rbind(res.pos,res.neg)
  
  gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  res <- merge(res,gmm.names,by.x="pathway", by.y="Module", all.x=T)
  res[,exposure:=toupper(exposure)]
  
  
  res <- res[,.(correlation,exposure,pathway,Name,HL1,HL2,pval,q.value,ES,NES,size)]
  
  setnames(res,c("pathway","pval","q.value"),c("Gut metabolic module","p.value","FDR_p"))
  
  setDF(res)
  
  # wb = loadWorkbook("Supp.tables.xlsx")
  # removeSheet(wb, sheetName = "GMM enrichment analysis")
  # yourSheet <- createSheet(wb, sheetName="GMM enrichment analysis")
  # addDataFrame(res, yourSheet, col.names=T,row.names=F)
  # saveWorkbook(wb, file = "Supp.tables.xlsx")
 
  write.xlsx2(res, "Supp.tables.xlsx", sheetName="GMM enrichment analysis", col.names=T,
              row.names=F, append=T)
  
  # Suppl Table 9 - GMM enrichment analysis - sensitivity analysis ####
  
  message("GMM enrichment - SA")
  
  res.pos <- fread(paste0(input,"ea_GMM_pos.tsv"))
  cols <- c("NES","pval","q.value")
  res.pos[,(cols):=lapply(.SD,round.large),.SDcols=cols]
  
  bmi.pw <- res.pos[exposure=="BMI" & q.value<.05,pathway]
  
  res.pos_sa <- fread("/home/baldanzi/Sleep_apnea/Results/ea_GMM_pos_sa_metformin.tsv")
  res.pos_sa[,exposure:=toupper(exposure)]
  res.pos_sa[,(cols):=lapply(.SD,round.large),.SDcols=cols]
  
  ahi.pw <- res.pos[exposure=="ahi" & q.value<.05 & !pathway %in% bmi.pw,
                    .(pathway,pval,q.value,NES)]
  
  names(ahi.pw) <- paste0(names(ahi.pw), "_main")
  
  ahi.pw_sa <- res.pos_sa[exposure=="AHI" & pathway %in% ahi.pw$pathway,
                    .(pathway,exposure,pval,q.value,NES)]
  
  ahi.pw_sa <- merge(ahi.pw_sa, ahi.pw, by.x="pathway",by.y="pathway_main", all.x=T)
  
  ahi.pw_sa <- merge(ahi.pw_sa, gmm.names, by.x="pathway",by.y="Module", all.x=T)
 
  names(ahi.pw_sa) <- gsub("q.value","FDR_p", names(ahi.pw_sa))
  names(ahi.pw_sa) <- gsub("pval","p.value", names(ahi.pw_sa))
  
  setcolorder(ahi.pw_sa,c("exposure","pathway",'Name','HL1','HL2',"p.value_main","FDR_p_main","NES_main"))
  
  
  
  t90.pw <- res.pos[exposure=="t90" & q.value<.05 & !pathway %in% bmi.pw,
                    .(pathway,pval,q.value,NES)]
  
  names(t90.pw) <- paste0(names(t90.pw), "_main")
  
  t90.pw_sa <- res.pos_sa[exposure=="T90" & pathway %in% t90.pw$pathway,
                          .(pathway,exposure,pval,q.value,NES)]
  
  t90.pw_sa <- merge(t90.pw_sa, t90.pw, by.x="pathway",by.y="pathway_main", all.x=T)
  
  t90.pw_sa <- merge(t90.pw_sa, gmm.names, by.x="pathway",by.y="Module", all.x=T)
  
  
  names(t90.pw_sa) <- gsub("q.value","FDR_p", names(t90.pw_sa))
  names(t90.pw_sa) <- gsub("pval","p.value", names(t90.pw_sa))
  
  setcolorder(t90.pw_sa,c("exposure","pathway","Name", "HL1", "HL2",
                          "p.value_main","FDR_p_main","NES_main"))
  
  
  break1 <- data.table(exposure="exposure",
                       pathway="pathway",
                       Name="",
                       HL1="",
                       HL2="",
                       p.value_main = "Main analysis",
                       FDR_p_main = "Main analysis",
                       NES_main = "Main analysis",
                       p.value = "Removed metformin users",
                       FDR_p = "Removed metformin users",
                       NES = "Removed metformin users")
  
  

  res <- rbind(break1,ahi.pw_sa,break1,t90.pw_sa)
  
   # wb = loadWorkbook("Supp.tables.xlsx")
   # removeSheet(wb, sheetName = "GMM enrichment_sens. analysis")
   # yourSheet <- createSheet(wb, sheetName = "GMM enrichment_sens. analysis")
   # addDataFrame(res, yourSheet, col.names=T,row.names=F)
   # saveWorkbook(wb, file = "Supp.tables.xlsx")
  
  write.xlsx2(res, "Supp.tables.xlsx", sheetName="GMM enrichment_sens. analysis", col.names=T,
              row.names=F, append=T)
  
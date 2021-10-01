# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-09-24

# Script to describe and compare self-reported sleep apnea with
# sleep apnea diagnosed with ApneaLink

# Analysis is restricted to Uppsala to avoid comparison between sites and focus
# on comparison between the sleep apnea measured vs self-reported 

rm(list=ls())

library(data.table)
library(compareGroups)
library(sjmisc)
library(vegan)


  pheno=fread("/home/baldanzi/Datasets/sleep_SCAPIS/SCAPIS-DATA-PETITION-170-20210315.csv", header = T, na.strings=c("", "NA"))
  pheno = pheno[Site=="Site5", ]
  setnames(pheno, "Subject", "SCAPISid") #renaming the id variable

  data.MGS = fread("/home/baldanzi/Datasets/MGS/clean/MGS_relative_abundance_4839_upp_4980_malmo.tsv")
  data.MGS[SCAPISid=="", SCAPISid:=id]
  data.MGS = data.MGS[SCAPISid %in% grep("5-",data.MGS$SCAPISid, value = T ),] 
  
  #Merge
  pheno <- merge(pheno, data.MGS, by="SCAPISid", all.y=T)
  pheno[,self.apnea:=cqhe058]
  pheno[self.apnea=="NO_ANSWER", self.apnea:=NA]


  #Age
  pheno[,age:=AgeAtVisitOne]
  
  #Sex
  pheno$Sex = factor(pheno$Sex, levels = c("MALE", "FEMALE"), 
                     labels = c("male", "female"))
  
  #smoking 
  setnames(pheno, 'derived_smoke_status', 'smokestatus')
  pheno$smokestatus = rec(pheno$smokestatus, 
                          rec = "NEVER=never; CURRENT=current; EX_SMOKER=former; NA,UNKNOWN=NA")
  pheno$smokestatus = factor(pheno$smokestatus, levels = c("never", "former", "current"))
  
  # Highest education 
  pheno$educat = rec(pheno$cqed001, rec = "-99=NA; else=copy")
  pheno$educat = factor(pheno$educat, levels = c(0,1,2,3), 
                        labels = c("uncompleted primary or lower secondary", 
                                   "lower secondary education", "upper secondary education", 
                                   "university education"))
  
  # Self-reported physical activity 
  pheno$leisurePA = rec(pheno$cqpa012, rec = "4=NA; else=copy")
  pheno$leisurePA = factor(pheno$leisurePA, levels = c(0,1,2,3), 
                           labels = c("mostly sedentary", "moderate activity", 
                                      "regular and moderate activity ", "regular exercise or training"))
  
  # Diabetes variable
  pheno$diabd = rec(pheno$Diabetes, rec = "NORMOGLYCEMIA=0; ELEV_HBA1C,IFG=1; KNOWN_DM,NEW_DM=2")
  pheno$diabd = factor(pheno$diabd, levels = c(0,1,2), 
                       labels = c("normoglycemic", "impaired glucose tolerance", "type 2 diabetes" ))
  # Self-reported hypertension
  pheno$hypertension = factor(pheno$cqhe034, levels = c("NO", "YES"), labels = c("no", "yes"))
  
  # Self-reported dyslipidemia 
  pheno$dyslipidemia = factor(pheno$cqhe036, levels = c("NO", "YES"), labels = c("no", "yes"))
  
  
  # Self-reported diabetes medication 
  pheno$diabmed = factor(pheno$cqhe039, levels = c("NO", "YES"), labels = c("no", "yes"))
  
  # Self-reported hypertension medication 
  pheno$hypermed = factor(pheno$cqhe035,  levels = c("NO", "YES"), labels = c("no", "yes"))
  
  # Self-reported dyslipidemia medication 
  pheno$dyslipmed = factor(pheno$cqhe037,  levels = c("NO", "YES"), labels = c("no", "yes"))
  
  # Self-reported Sleep apnea treatment  
  pheno$apneatto_self = factor(pheno$cqhe060,  levels = c("NO", "YES"), labels = c("no", "yes"))
  
  # Self-reported Sleep apnea treatment CPAP
  pheno$cpap_self = ifelse(pheno$cqhe061=="CPAP", "yes", "no")
  
  
  
  # Import ApneaLink data ####
  sleep <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording/sleep.rds")
  
  # Merge
  pheno <- merge(pheno, sleep, by='SCAPISid', all.x=T, all.y=F)
  
  
  # Shannon 
  mgs.names = grep("____",names(pheno),value=T) # vector with MGS names 
  pheno[,shannon:=diversity(pheno[,mgs.names, with=F],index="shannon")]
  
  
  # Comparison table 
  
  pheno[,type.apnea:=NULL]
  pheno[self.apnea=='YES',type.apnea:=1]
  pheno[self.apnea=='YES' & cpap_self!="yes",type.apnea:=2]
  pheno[OSAcat=="Moderate",type.apnea:=3]
  pheno[OSAcat=="Severe",type.apnea:=4]
  
  pheno[,type.apnea:= factor(type.apnea, levels = 1:4, 
                             labels=c("Self-rep.", "Self-rep.No tt",
                                      "Moderate", "Severe"))]
  
    t <- compareGroups(type.apnea~age+Sex+BMI+smokestatus+Alkohol+educat+leisurePA+
                         diabd+hypertension+dyslipidemia+diabmed+hypermed+dyslipmed+
                         apneatto_self+cpap_self+ahi+t90+shannon, data=pheno)
    t1 <- createTable(t)
    
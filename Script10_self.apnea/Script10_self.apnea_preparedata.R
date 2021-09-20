# Sleep apnea and microbiota project 

# Gabriel Baldanzi - 2021-09-17

# Last update: 2021-09-17

# Script to prepare data for the Self-apnea analysis 

  library(data.table)
  library(sjmisc)
  library(vegan)

  rm(list=ls())

all.pheno=fread("/home/baldanzi/Datasets/sleep_SCAPIS/SCAPIS-DATA-PETITION-170-20210315.csv", header = T, na.strings=c("", "NA"))
setnames(all.pheno, "Subject", "SCAPISid") #renaming the id variable

data.MGS = fread("/home/baldanzi/Datasets/MGS/clean/MGS_relative_abundance_4839_upp_4980_malmo.tsv")
data.MGS[SCAPISid=="", SCAPISid:=id]

all.pheno <- merge(all.pheno, data.MGS, by="SCAPISid", all.y=T)

all.pheno[,self.apnea:=cqhe058]
all.pheno[self.apnea=="NO_ANSWER", self.apnea:=NA]



#Age
all.pheno[,age:=AgeAtVisitOne]

#Sex
all.pheno[,Sex:= factor(Sex, levels = c("MALE", "FEMALE"), 
                          labels = c("male", "female"))]

# Smoking
setnames(all.pheno, 'derived_smoke_status', 'smokestatus')
all.pheno$smokestatus = rec(all.pheno$smokestatus, 
                        rec = "NEVER=never; CURRENT=current; EX_SMOKER=former; NA,UNKNOWN=NA")
all.pheno$smokestatus = factor(all.pheno$smokestatus, levels = c("never", "former", "current"))


# Clinical microbiomics variables 
cmvar = fread("/home/baldanzi/Datasets/MGS/CM_lab_variables_4838_upp_4980_malmo.tsv",header=T, na.strings=c("", "NA")) 
selected_id = fread("/home/baldanzi/Datasets/scapis_idkey/selected_participants_id_sample_id.tsv", 
                    header=T, na.strings=c("", "NA"))
cmvar = merge(selected_id, cmvar, by = "sample.id", all=F, all.x=T )
ids=fread("/home/baldanzi/Datasets/scapis_idkey/id_conversion.txt", header=T) 
ids$id = as.character(ids$export_id)
cmvar = merge(cmvar, ids, by = "id", all=F, all.x=T)
cmvar[is.na(subject_id), subject_id:=id]
setnames(cmvar, "subject_id", "SCAPISid")

all.pheno <- merge(all.pheno, cmvar, by="SCAPISid", all.x=T)

  # Calculate shannon index 
  a = grep("____",names(all.pheno),value=T) # vector with MGS names 
  all.pheno[,shannon:=diversity(all.pheno[,a, with=F],index="shannon")]
  
  
  # Remove individuals with sleep records 
  only.upp <- readRDS(file="/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

  all.pheno <- all.pheno[-which(SCAPISid %in% only.upp[valid.t90=='yes',]),SCAPISid]
  

  saveRDS(all.pheno, file='/home/baldanzi/Datasets/sleep_SCAPIS/self.apnea.rds')
  
  apnea.indiv <- all.pheno[self.apnea=="YES",]
  apnea.indiv[,CPAP:='NO']
  apnea.indiv[cqhe061=='CPAP',CPAP:='YES']
  
  

  
  
  saveRDS(apnea.indiv, file='/home/baldanzi/Datasets/sleep_SCAPIS/apnea.indiv.rds')




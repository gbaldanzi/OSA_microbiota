# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will calculate alpha diversity in those individuals that have valid.ahi 

rm(list=ls())


# Uploading dataset with variables of class "factor"
load("/home/baldanzi/Datasets/sleep_SCAPIS/data_table1")
# Importing data on individuals with valid.ahi measurement 
valid.ahi = fread("/home/baldanzi/Datasets/sleep_SCAPIS/validsleep.MGS.Upp.tsv", header=T, na.strings=c("", "NA"))
#merging with the data set that has variables as factor 
a = names(dat1[,-"SCAPISid"])
valid.ahi[,(a):=NULL]
valid.ahi = merge(dat1,valid.ahi, by = "SCAPISid", all=T, all.x=F, all.y=F)

valid.ahi[,visit.month:=as.factor(valid.ahi$visit.month)]

#Calculate Shannon diversity ####
  a = grep("____",names(valid.ahi),value=T) # vector with MGS names 
  valid.ahi$shannon=diversity(valid.ahi[,a, with=F],index="shannon") 

  fwrite(valid.ahi, file = "/home/baldanzi/Datasets/sleep_SCAPIS/validsleep_MGS.shannon_Upp.tsv", sep = "\t")
  
  saveRDS(valid.ahi, file = "/home/baldanzi/Datasets/sleep_SCAPIS/validsleep_MGS.shannon_Upp.rds")

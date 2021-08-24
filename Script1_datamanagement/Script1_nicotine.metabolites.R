# Project: Sleep apnea and gut microbiota

# Nicotine could be considered as an alternative variable for smoking status

# Here is a dataset to extract information on nicotine metabolites and merge with 
# the smoking status variable 

library(data.table)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
# Importing data on individuals with valid.ahi measurement (phenotype+MGS)
valid.ahi = fread("validsleep.MGS.Upp.tsv", header=T, na.strings=c("", "NA"))

# Nicotine metabolites ("cotinine", "hydroxycotinine", "cotinine N-oxide") ####
metabolon=fread("/home/baldanzi/Datasets/Metabolon/clean/scapis_metabolon_batchnorm_scapisid.tsv", header=T)

a = c("SCAPISid" ,"MET_848","MET_100002717", "MET_100002719")
nicotine=metabolon[,a,with=F]

# Merging nicotine metabolites with smokestatus 
a = c("SCAPISid", "smokestatus", "cqto015")
nicotine = merge(nicotine, valid.ahi[,a, with=F], by="SCAPISid", all.x = F, all.y=T)
setnames(nicotine,"cqto015", "snus_ever1mo")

# Saving nicotine metabolites data 
fwrite(nicotine, file = "/home/baldanzi/Datasets/Metabolon/nicotine_metab.csv", sep=",")
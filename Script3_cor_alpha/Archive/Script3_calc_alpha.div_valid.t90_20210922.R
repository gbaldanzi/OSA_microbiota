# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will calculate alpha diversity in those individuals that have valid 
# T90 measurment 

rm(list=ls())

# Loading packages
require(vegan)
require(data.table)

# Importing data
valid.t90 = fread("/home/baldanzi/Datasets/sleep_SCAPIS/validodi.MGS.Upp.tsv", sep = "\t")
valid.t90 = valid.t90[!is.na(sat90),]
setnames(valid.t90, "sat90", "t90")

# Create factor variables 
valid.t90[,visit.month:=as.factor(valid.t90$visit.month)]
valid.t90[,plate:=as.factor(valid.t90$plate)]

# Create a grouping variable of t90
valid.t90[t90<10,t90cat :=1]
valid.t90[t90>=10 & t90<20,t90cat :=2]
valid.t90[t90>=20,t90cat :=3]
valid.t90[,t90cat:= factor(valid.t90$t90cat, levels = c(1,2,3),
                             labels = c("<10","10-20",">20"))]

#Calculate Shannon diversity ####
a = grep("____",names(valid.t90),value=T) # vector with MGS names 
valid.t90$shannon=diversity(valid.t90[,a, with=F],index="shannon") #estimating shannon per individual


fwrite(valid.t90, file = "/home/baldanzi/Datasets/sleep_SCAPIS/valid.t90_MGS.shannon_Upp.tsv", sep = "\t")

saveRDS(valid.t90, file = "/home/baldanzi/Datasets/sleep_SCAPIS/valid.t90_MGS.shannon_Upp.rds")
#-----------------------------------------------------------------------------#

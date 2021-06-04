# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics continued : script to DESCRIBE alpha and beta-diversity

# This script will calculate alpha diversity and create plots and table describing alpha
#diversity in relation AHI and OSA, as well as other covarites. 

# Descriptive Statistics 

rm(list=ls())

# Loading packages
library(tidyverse)
library(data.table)
library(vegan)
library(Hmisc)
library(compareGroups)
library(flextable)
library(pheatmap)
library(RColorBrewer)
library(robCompositions)
library(grid)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")

# Uploading dataset with variables of class "factor"
load("data_table1")
# Importing data on individuals with valid.ahi measurement 
valid.ahi = fread("validsleep.MGS.Upp.tsv", header=T, na.strings=c("", "NA"))
#merging with the data set that has variables as factor 
a = names(dat1[,-"SCAPISid"])
valid.ahi[,(a):=NULL]
valid.ahi = merge(dat1,valid.ahi, by = "SCAPISid", all=T, all.x=F, all.y=F)


#Shannon diversity ####
a = grep("____",names(valid.ahi),value=T) # vector with MGS names 
valid.ahi$shannon=diversity(valid.ahi[,a, with=F],index="shannon") #estimating shannon per individual

fwrite(valid.ahi, file = "validsleep_MGS.shannon.BC_Upp.tsv", sep = "\t")



# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will calculate the aitchison distance (euclidean distance of clr 
# transformed counts) between the participants with valid t90 


rm(list=ls())

# Loading packages
library(ape)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")


# Import data
valid.t90 <- readRDS("valid.t90_MGS.shannon_Upp.rds")

# AITCHISON DISTANCES ####

# Importing count data:
  count = readRDS('/home/baldanzi/Datasets/MGS/clean/upugut03.mgsCounts_clean.rds')

# Zeros are replaced with 1 (there were no previous 1 in the data)
  count[count==0]=1

# Central log-transformation of count data 
  count.clr = cenLR(count)
  a = rownames(count.clr$x.clr)
  count = as.data.frame(count.clr$x.clr)
  rownames(count) = a
  count$SCAPISid = a

# Only keeping the individuals that have valid T90
  count=count[count$SCAPISid %in% valid.t90[,SCAPISid],]

# Merging counts with valid.t90 
  a = c("SCAPISid", "t90cat", "ahi", "age", "Sex")
  count=merge(count,valid.t90[,a,with=F], by="SCAPISid",all=F, all.x=F, all.y=T)

# Aitchison distance ####
  mgs_count = as.matrix(count[,grep("HG3A",names(count))])
  rownames(mgs_count) <- count$SCAPISid

  t0 = Sys.time()
  atdist=vegdist(mgs_count, method="euclidean")
  t1 = Sys.time()
  print(paste0("Vegan time=" ,t1-t0))


# Transform aitchison distances into a matrix
atdist_matrix=as.matrix(atdist)

# Save Aitchison distance matrix  
fwrite(atdist_matrix,
       file = "/home/baldanzi/Datasets/sleep_SCAPIS/T90.ADmatrix.csv",sep=',')

# PCoA of Aitchison distances ####
pcoa.ad <- pcoa(atdist_matrix)
save(pcoa.ad, file = 'pc_AD_t90')


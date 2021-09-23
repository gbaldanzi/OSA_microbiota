# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will calculate the aitchison distance (euclidean distance of clr 
# transformed counts) between the participants with valid ahi 


rm(list=ls())

# Loading packages
library(ape)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")

# Import pheno data
valid.ahi <- readRDS("validsleep_MGS.shannon_Upp.rds")

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
  
# Only keeping the individuals that have valid sleep recording
  count=count[count$SCAPISid %in% valid.ahi[,SCAPISid],]

# Merging counts with valid.ahi 
  a = c("SCAPISid", "OSAcat", "ahi", "age", "Sex")
  count=merge(count,valid.ahi[,a,with=F], by="SCAPISid",all=F, all.x=F, all.y=T)

# Aitchison distance ####
  rownames(count)=count$SCAPISid
  mgs_count = as.matrix(count[,grep("HG3A",names(count))])


 t0 = Sys.time()
 atdist=vegdist(mgs_count, method="euclidean")
 t1 = Sys.time()
 print(paste0("Vegan time=" ,t1-t0))

 
# Transform aitchison distances into a matrix
 atdist_matrix=as.matrix(atdist)
 
# Save Aitchison distance matrix  
 fwrite(atdist_matrix,
        file = "/home/baldanzi/Datasets/sleep_SCAPIS/OSA.aitchison_distmatrix.csv",sep=',')
 
 # PCoA of Aitchison distances ####
 pcoa.ad <- pcoa(atdist_matrix)
 save(pcoa.ad, file = 'pc_AD')

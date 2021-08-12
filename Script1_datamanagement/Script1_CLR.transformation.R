# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script to produce a data table of central log-ratio transformation to count data 

rm(list=ls())

# Loading packages
library(ape)

  # Import count data 
  count = readRDS('/home/baldanzi/Datasets/MGS/clean/upugut03.mgsCounts_clean.rds')
  
  # Zeros are replaced with 1 (there were no previous 1 in the data)
  count[count==0]=1
  
  # Central log-transformation of count data 
  count.clr = cenLR(count)
  a = rownames(count.clr$x.clr)
  count = as.data.frame(count.clr$x.clr)
  rownames(count) = a
  
  
  # Adjust MGS names ( from HG3A to names____HG3A)
  input.wrf <-  "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"
  taxonomy <-  fread(paste0(input.wrf,"taxonomy/MGS_taxonomic_information.tsv"))
  taxonomy <-  taxonomy[,.(maintax_mgs, mgs)]
  
  newnames <- taxonomy[match(names(count),mgs), maintax_mgs]
  
  names(count) <- newnames
  
  count$SCAPISid = a
  
  setDT(count)
  setcolorder(count, neworder = "SCAPISid")
  
  fwrite(count,file = '/home/baldanzi/Datasets/MGS/clean/MGS_clr.transformed_4839_upp.tsv')
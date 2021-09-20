# Script 3 - source

# Run this script to calculate alpha and beta diversity as well as create 
# the respective plots 

library(tidyverse)
library(data.table)
library(vegan)
library(robCompositions)
library(pheatmap)

#1
print("start - script3_calc_alpha.div_valid.ahi")
setwd('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
  source("Script3_calc_alpha.div_valid.ahi.R")
#2  
print("start - Script3_calc_alpha.div_valid.ahi_plots.R")
setwd('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
  source("Script3_calc_alpha.div_valid.ahi_plots.R")

#3
print("start - Script3_calc_alpha.div_valid.t90.R")
  setwd('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
  source("Script3_calc_alpha.div_valid.t90.R")
#4
print("start - Script3_calc_alpha.div_valid.t90_plots.R")
setwd('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
  source("Script3_calc_alpha.div_valid.t90_plots.R")

#5
print("start - Script3_calc_mgs.prevalence.R")  
setwd('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
#  source("Script3_calc_mgs.prevalence.R")

#6 
message("start - Script3_calc_alpha.div_pheno.R")
setwd('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
  source("Script3_calc_alpha.div_pheno.R")

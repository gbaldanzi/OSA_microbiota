# Script 3 - source

# Run this script to calculate Bray-curtis dissimilary matrix and produce PCoA plots

  library(tidyverse)
  library(data.table)
  library(vegan)
  library(ape)
  library(RColorBrewer)
  library(grid)
  library(pheatmap)

  input <- '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/Script3_calc_beta/'

  # Import data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

#1
  message("start - Script3_calc_bc_ahi.R")  
  source(paste0(input,"Script3_calc_bc_ahi.R")) # worked
  
#2
  message("start - Script3_calc_bc_ahi_plots.R") 
  source(paste0(input,"Script3_calc_bc_ahi_plots.R")) #worked

#3
  #message("start - Script3_calc_bc_t90.R")
  #source(paste0(input,"Script3_calc_bc_t90.R")) #worked 
  
#4 
  #message("start - Script3_calc_bc_t90_plots.R")
  #source(paste0(input,"Script3_calc_bc_t90_plots.R"))
  

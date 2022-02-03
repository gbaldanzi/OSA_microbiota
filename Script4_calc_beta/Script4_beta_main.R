# Script 3 - source

# Run this script to calculate Bray-curtis dissimilary matrix and produce PCoA plots

  library(tidyverse)
  library(data.table)
  library(vegan)
  library(ape)
  library(RColorBrewer)
  library(grid)
  
  # Folder where scripts are located 
    input.script <- 'Script4_calc_beta/'
    
  # Folder where dataset is located 
    input <- "/home/baldanzi/Datasets/sleep_SCAPIS/"
    
  # Folder for output results 
    output.plot <- "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Import data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

  # Calculate Bray-curtis for participants with valid AHI 
   
  source(paste0(input.script,"Script4_calc_bc_valid.ahi.R")) 
  
  # Calculate Bray-curtis for participants with valid T90 and ODI
 
  source(paste0(input.script,"Script4_calc_bc_valid.t90.R")) 
  
  # PCoA for AHI, T90, and ODI severity groups 
  
  source(paste0(input.script, "Script4_merge.plots.R"))
  

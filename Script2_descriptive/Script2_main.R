# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# Last update: 2022-01-02

# Loading packages
library(tidyverse)
library(data.table)
library(Hmisc)
library(compareGroups)


  # Import data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  
  #message("Run Descriptive - Sleep records")
  #source('Script2_descriptive/Script2_descriptive_sleeprecords.R')
  
  message("Run Descriptive - Phenotypes")
  source('Script2_descriptive/Script2_descriptivestatistics.R')
  

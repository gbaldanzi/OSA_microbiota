# Data management update 

# This script is used to fast update the datasets after changes made 
# on the data management script 

library(tidyverse)
library(data.table)
library(vegan)

source('Script1_sleeprecords.R')

setwd('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
source('Script1_datamanagement.R')

setwd('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
source("Script3_calc_alpha.div_valid.ahi.R")

setwd('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
source("Script3_calc_alpha.div_valid.t90.R")

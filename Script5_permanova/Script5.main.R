# Script 5 - Main 

# Script to run the permutation analysis between OSA/AHI/T90 and BC/AD

# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)
 

  source('perma_osa_bc/perma_osa_bc_main.R')

  source('perma_ahi_bc/perma_ahi_bc_main.R')

  source('perma_t90_bc/perma_t90_bc_main.R')

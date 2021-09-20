# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-31

# Scatter plot of MGS correlations coefficients for AHI and BMI, as well as 
# for T90 and BMI. 

  pacman::p_load(data.table, tidyverse)

  # Import data
  res = fread(paste0(output, "cor_all.var_mgs.tsv"))
  
  res$q.value[res$q.value>=0.001] <- round(res$q.value[res$q.value>=0.001] , digits = 3)
  
  mgs.t90 <- res[exposure=="t90" & q.value<0.05, MGS]
  mgs.ahi <- res[exposure=="ahi" & q.value<0.05, MGS]
  mgs.bmi <- res[exposure=="BMI" & q.value<0.05, MGS]
  mgs.ahionly <- mgs.ahi[!mgs.ahi %in% mgs.bmi]
  mgs.t90only <- mgs.t90[!mgs.t90 %in% mgs.bmi]
  

  # Prepare plot data
  data.plot <- res %>% select("MGS", "exposure", "cor.coefficient") %>% 
    pivot_wider(id_cols = MGS, names_from = exposure, values_from = cor.coefficient)
  
  data.plot <- data.plot[data.plot$MGS %in% mgs.m1,]
  
  data.plot$ahionly <- 0
  data.plot$ahionly[data.plot$MGS %in% mgs.ahionly] <- 1
  
  data.plot$t90only <- 0
  data.plot$t90only[data.plot$MGS %in% mgs.t90only] <- 1
  
  p <- ggplot(data.plot,aes(x=ahi, y=BMI,color = ahionly)) +
    geom_point()
  p
  ggsave('cor.coef_ahi_bmi.jpeg')
  
  p <- ggplot(data.plot,aes(x=t90, y=BMI,color = t90only)) +
    geom_point()
  p
  ggsave('cor.coef_t90_bmi.jpeg')
  
  
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-01-30

# Inferential Statistics 

# This code will produce a Venn Diagram from the correlations of AHI, BMI, 
# and T90 with MGSs from the Sensitivity analysis 

  rm(list = ls())
  # Loading packages 
  pacman::p_load(data.table,ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results 
  res.ahi <- fread(paste0(input,"corsa_ahi_mgs.tsv"))
  res.bmi <- fread(paste0(input,"corsa_BMI_mgs.tsv"))
  res.t90 <- fread(paste0(input,"corsa_t90_mgs.tsv"))
  
  res.list = list(res.ahi,res.bmi,res.t90)
  
  
  # filter MGS significant at the FDR p-value<0.05
  res.list <- lapply(res.list, function(x){x[q.value>=0.001, q.value:=round(q.value, digits = 3)]})
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})
  
  names(mgs.fdr) <- c("AHI","BMI","T90")
  
  # Create VennDiagram    
  venn.m2 <- ggvenn(mgs.fdr,
                    fill_color = c("yellow","orange","pink"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .5) +
    ggtitle("MGS correlated with AHI, T90% and BMI - sens anlaysis") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_bmiahi_sa.png",plot = venn.m2,device = "png", path=output.plot)
  





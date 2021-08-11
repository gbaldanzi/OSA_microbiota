# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-18

# Inferential Statistics 

# This code will produce a Venn Diagram from the correlations of AHI, BMI, 
# and T90 with MGSs in model 1

  rm(list = ls())
  # Loading packages 
  pacman::p_load(data.table,ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results 
  res.ahi <- fread(paste0(input,"cor_ahi_mgs.tsv"))
  res.bmi <- fread(paste0(input,"cor_BMI_mgs.tsv"))
  res.t90 <- fread(paste0(input,"cor_t90_mgs.tsv"))
  
  res.list = list(res.ahi,res.bmi,res.t90)
  
  
  # filter MGS significant at the FDR p-value<0.05
  res.list <- lapply(res.list, function(x){x[q.value>=0.001, q.value:=round(q.value, digits = 3)]})
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})
  
  names(mgs.fdr) <- c("AHI","BMI","T90")
  
  # Create VennDiagram    
  venn.m1 <- ggvenn(mgs.fdr,
                    fill_color = c("yellow","orange","pink"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .5) +
    ggtitle("MGS correlated with AHI, T90% and BMI - model 1") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_bmiahi_m1.png",plot = venn.m1,device = "png", path=output.plot)
  
  # All MGS model2 
  mgs.m1 <-  c(mgs.fdr[[1]], mgs.fdr[[2]], mgs.fdr[[3]])
  mgs.m1 <- unique(mgs.m1)
  
  # Save MGS that were related to either AHI, T90 or BMI. 
  # These MGS will be used to the second model. 
  saveRDS(mgs.m1, file = '/home/baldanzi/Sleep_apnea/Results/mgs.m1.rds')




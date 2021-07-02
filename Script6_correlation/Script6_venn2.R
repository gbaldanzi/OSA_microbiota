# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-18

# Inferential Statistics 

# This code will produce a Venn Diagram from the correlations of AHI, BMI, 
# and T90 with MGSs in model 2

  rm(list = ls())
  # Loading packages 
  pacman::p_load(data.table,ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results 
  res.ahi <- fread(paste0(input,"cor2_ahi_mgs.tsv"))
  res.bmi <- fread(paste0(input,"cor2_BMI_mgs.tsv"))
  res.t90 <- fread(paste0(input,"cor2_t90_mgs.tsv"))
  
  res.list = list(res.ahi,res.bmi,res.t90)
  
  
  # filter MGS significant at the FDR p-value<0.05
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})
  
  names(mgs.fdr) <- c("AHI","BMI","T90")
  
  # Create VennDiagram    
  venn.m2 <- ggvenn(mgs.fdr,
                    fill_color = c("yellow","orange","pink"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .5) +
    ggtitle("MGS correlated with AHI, T90% and BMI - model 2") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_bmiahi_m2.png",plot = venn.m2,device = "png", path=output.plot)
  
  # All MGS model1 
  mgs.m2 <-  c(mgs.fdr[[1]], mgs.fdr[[2]], mgs.fdr[[3]])
  mgs.m2 <- unique(mgs.m2)
  
  # Save MGS that were related to either AHI, T90 or BMI. 
  # These MGS will be used to the second model. 
  saveRDS(mgs.m2, file = '/home/baldanzi/Sleep_apnea/Results/mgs.m2.rds')




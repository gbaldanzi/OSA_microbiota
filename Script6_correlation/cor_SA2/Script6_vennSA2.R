# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-09-02

# Inferential Statistics 

# This code will produce a Venn Diagram from the correlations of AHI, BMI, 
# and T90 with MGSs in sensitivity analysis 2

  rm(list = ls())
  # Loading packages 
  pacman::p_load(data.table,ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results 
  res <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  
  res.bmi <-  res[exposure=="BMI",]
  res.ahi <- fread(paste0(input,"corSA2_ahi_mgs.tsv"))
  res.t90 <- fread(paste0(input,"corSA2_t90_mgs.tsv"))
  
  res.list <- list(AHI = res.ahi, 
                   T90 = res.t90, 
                   BMI= res.bmi)

  # filter MGS significant at the FDR p-value<0.05
  res.list <- lapply(res.list, function(x){x[q.value>=0.001, q.value:=round(q.value, digits = 3)]})
  
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})
  
   # Create VennDiagram    
  venn.sa2 <- ggvenn(mgs.fdr,
                    fill_color = c("yellow","orange","pink"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .5) +
    ggtitle("MGS correlated with AHI, T90%, and BMI - SA2") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  
    ggsave("Venn_bmiahi_SA2.png",plot = venn.sa2,device = "png", path=output.plot)
  

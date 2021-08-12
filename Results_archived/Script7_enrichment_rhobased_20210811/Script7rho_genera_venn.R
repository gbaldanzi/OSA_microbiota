# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-05

# Inferential Statistics 

# This code will produce a Venn Diagram from the overrepresentation analysis of genera
# correlated to AHI, T90 and BMI

# Rho based 

  rm(list = ls())
  # Loading packages 
  pacman::p_load(data.table,ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results
  res.both <- readRDS(file = paste0(input,"ea_genera_rho.rds"))

  # genera list 
  genera = lapply(res.both,function(x){x[q.value<0.05,genera]})
  
  names(genera) <- c("AHI","BMI","T90")
  
  # Create VennDiagram    
  venn1 <- ggvenn(genera,
                  fill_color = c("yellow","orange","pink"), 
                  stroke_color = "white",
                  stroke_size = .2, 
                  show_percentage = F, 
                  fill_alpha = .5) +
    ggtitle("Genera overrepresented in the AHI, T90, and BMI both correlations") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_genera_rho.png",plot = venn1,device = "png", path=output.plot)
  

  

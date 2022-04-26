# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-24

# Inferential Statistics 

# This code will produce a Venn Diagram from the overrepresentation analysis of genera
# correlated to AHI, T90 and BMI

  rm(list = ls())
  # Loading packages 
  pacman::p_load(data.table,ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results
  res.both <- readRDS(file = paste0(input,"ea_genera_both.rds"))
  res.pos <- readRDS(file = paste0(input,"ea_genera_pos.rds"))
  res.neg <- readRDS(file = paste0(input,"ea_genera_neg.rds"))
  
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
  ggsave("Venn_genera_both.png",plot = venn1,device = "png", path=output.plot)
  

  # POSITIVE correlation 
  # genera list 
  genera = lapply(res.pos,function(x){x[q.value<0.05,genera]})
  
  names(genera) <- c("AHI","BMI","T90")
  
  # Create VennDiagram    
  venn2 <- ggvenn(genera,
                    fill_color = c("yellow","orange","pink"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .5) +
    ggtitle("Genera overrepresented in the AHI, T90, and BMI positive correlations") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_genera_pos.png",plot = venn2,device = "png", path=output.plot)
  
  # NEGATIVE CORRELATIONS
  # genera list 
  genera = lapply(res.neg,function(x){x[q.value<0.05,genera]})
  
  names(genera) <- c("AHI","BMI","T90")
  
  # Create VennDiagram    
  venn3 <- ggvenn(genera,
                    fill_color = c("yellow","orange","pink"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .5) +
    ggtitle("Genera overrepresented in the AHI, T90, and BMI negative correlations") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_genera_neg.png",plot = venn3,device = "png", path=output.plot)
  





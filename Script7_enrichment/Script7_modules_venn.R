# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-24

# Inferential Statistics 

# This code will produce a Venn Diagram from the overrepresentation analysis of KEGG modules
# correlated to AHI, T90 and BMI

  rm(list = ls())
  # Loading packages 
  pacman::p_load(data.table,ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results
  res.both <- fread(paste0(input,"ea_modules_both.tsv"))
  res.pos <- fread(paste0(input,"ea_modules_pos.tsv"))
  res.neg <- fread(paste0(input,"ea_modules_neg.tsv"))

  # function for modules list 
  neat.module.list.fun <- function(data){
    module.list = list(
                  data[exposure=="ahi" & q.value<0.05,pathway],
                  data[exposure=="BMI" & q.value<0.05,pathway],
                  data[exposure=="t90" & q.value<0.05,pathway]
          )
    names(module.list) <- c("AHI","BMI","T90")
    return(module.list)
  }
  
  # listing 
  list.res.both <- neat.module.list.fun(res.both)
  list.res.pos <- neat.module.list.fun(res.pos)
  list.res.neg <- neat.module.list.fun(res.neg)
  

  # Create VennDiagram    
  venn1 <- ggvenn(list.res.both,
                  fill_color = c("yellow","orange","pink"), 
                  stroke_color = "white",
                  stroke_size = .2, 
                  show_percentage = F, 
                  fill_alpha = .5) +
    ggtitle("Modules overrepresented in the AHI, T90, and BMI") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_modules_both.png",plot = venn1,device = "png", path=output.plot)
  

  # POSITIVE correlation 

  # Create VennDiagram    
  venn2 <- ggvenn(list.res.pos ,
                    fill_color = c("yellow","orange","pink"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .5) +
    ggtitle("Modules overrepresented in the AHI, T90, and BMI positive correlations") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_modules_pos.png",plot = venn2,device = "png", path=output.plot)
  
  # NEGATIVE CORRELATIONS
  
  # Create VennDiagram    
  venn3 <- ggvenn(list.res.neg ,
                    fill_color = c("yellow","orange","pink"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .5) +
    ggtitle("Modules overrepresented in the AHI, T90, and BMI negative correlations") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_modules_neg.png",plot = venn3,device = "png", path=output.plot)
  





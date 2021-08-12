# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-07-05

# Inferential Statistics 

# This code will produce a Venn Diagram from the overrepresentation analysis of subpathways 
# among metabolites correlated MGSs

  rm(list = ls())
  # Loading packages 
  pacman::p_load(data.table,ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results
  res <- fread(paste0(input,"ea_subpath.tsv"))


  # function for modules list 
  neat.module.list.fun <- function(data){
    module.list = list(
                  data[exposure=="AHI" & q.value<0.05,pathway],
                  data[exposure=="BMI" & q.value<0.05,pathway],
                  data[exposure=="T90" & q.value<0.05,pathway]
          )
    names(module.list) <- c("AHI","BMI","T90")
    return(module.list)
  }
  
  # listing 
  list.res <- neat.module.list.fun(res)


  # Create VennDiagram    
  venn1 <- ggvenn(list.res,
                  fill_color = c("yellow","orange","pink"), 
                  stroke_color = "white",
                  stroke_size = .2, 
                  show_percentage = F, 
                  fill_alpha = .5) +
    ggtitle("Subpathways overrepresented for the AHI, T90, and BMI") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_subpathways.png",plot = venn1,device = "png", path=output.plot)
  

# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-10-25

# Last update: 2021-10-25

# This code will produce a Venn Diagram from the overrepresentation analysis of GMM
# correlated to AHI, T90 and BMI

  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"
  #output.plot= "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512"
  
  
  # Importing results
  res.pos <- fread(paste0(input,"ea_GMM_pos.tsv"))
  res.neg <- fread(paste0(input,"ea_GMM_neg.tsv"))
  
  # function for modules list 
  neat.module.list.fun <- function(data){
    data[q.value>=0.001, q.value:=round(q.value, digits = 3)]
    module.list = list(
                  data[exposure=="ahi" & q.value<0.05,pathway],
                  data[exposure=="t90" & q.value<0.05,pathway],
                  data[exposure=="BMI" & q.value<0.05,pathway]
          )
    names(module.list) <- c("AHI","T90","BMI")
    return(module.list)
  }
  
  # listing 
  list.res.pos <- neat.module.list.fun(res.pos)
  list.res.neg <- neat.module.list.fun(res.neg)
  

  # POSITIVE correlation 

  # Create VennDiagram    
  venn2 <- ggvenn(list.res.pos ,
                    fill_color = c("orange","cornflowerblue","green4"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .6) +
    ggtitle("Pathways enriched among species \npositively correlated") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_GMM_pos.png",plot = venn2,device = "png", path=output.plot)
  
  # NEGATIVE CORRELATIONS
  
  # Create VennDiagram    
  venn3 <- ggvenn(list.res.neg ,
                    fill_color = c("orange","cornflowerblue","green4"), 
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .6) +
    ggtitle("GMM modules enriched among the MGSs negatively \ncorrelated to AHI, T90, or BMI") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_GMM_neg.png",plot = venn3,device = "png", path=output.plot)
  





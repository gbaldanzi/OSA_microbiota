# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-06-18

# This code will produce a Venn Diagram from the correlations of AHI, BMI, 
# and T90 with MGSs in model 2

  rm(list = ls())
  # Loading packages 
  pacman::p_load(data.table,ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results 
  res <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   BMI = res[exposure=="BMI",])

  # filter MGS significant at the FDR p-value<0.05
  res.list <- lapply(res.list, function(x){x[q.value>=0.001, q.value:=round(q.value, digits = 3)]})
  
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})
  
   # Create VennDiagram    
  venn.m2 <- ggvenn(mgs.fdr,
                    fill_color = c("orange","cornflowerblue","green4"),
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .6) +
    ggtitle("Microbiota species associated with AHI, T90%, and/or BMI\nfully adjusted model") +
    theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))
  
    ggsave("Venn_bmiahi_m2.png",plot = venn.m2,device = "png", 
           path=output.plot)
    ggsave("Venn_bmiahi_m2.pdf",plot = venn.m2,device = "pdf", 
           path='/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
    

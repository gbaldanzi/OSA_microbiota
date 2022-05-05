# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Venn diagram - Fig 2

# This code will produce a Venn Diagram from the FDR-significant correlations of AHI, ODI, 
# and T90 with species using the main model including and not including BMI

  rm(list = ls())
  # Loading packages 
  library(data.table)
  library(ggplot2)
  library(ggvenn)
  library(tidyr)
  library(cowplot)

  # results.folder and output folders 
  results.folder <- "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/"
  output.plot = "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/figures/"

  # Model not including BMI 
  res <- fread(paste0(results.folder,"cor_all.var_mgs.tsv"))
  res[q.value < 0.001, q.value := round(q.value,3)]

  res.list <- list(AHI = res[exposure=="ahi",],
                 T90 = res[exposure=="t90",],
                 ODI = res[exposure=="odi",])


  # filter MGS significant at the FDR p-value<0.05
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})

  # Create VennDiagram    
  venn.m1 <- ggvenn(mgs.fdr,
                  fill_color = c("orange","cornflowerblue","grey83"),
                  stroke_color = "white",
                  stroke_size = .2, 
                  show_percentage = F, 
                  fill_alpha = .6, 
                  set_name_size = 4,
                  text_size = 3) +
  ggtitle("Not adjusted for BMI") +
  theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))


  # Model including BMI ####

  # Importing results 
  res <- fread(paste0(results.folder,"cor.bmi_all.var_mgs.tsv"))
  res[q.value < 0.001, q.value := round(q.value,3)]

  res.list <- list(AHI = res[exposure=="ahi",],
                 T90 = res[exposure=="t90",],
                 ODI = res[exposure=="odi",])

  # filter MGS significant at the FDR p-value<0.05
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})

  # Create VennDiagram    
  venn.m2 <- ggvenn(mgs.fdr,
                  fill_color = c("orange","cornflowerblue","grey83"),
                  stroke_color = "white",
                  stroke_size = .2, 
                  show_percentage = F, 
                  fill_alpha = .6, 
                  set_name_size = 4,
                  text_size = 3) +
  ggtitle("BMI adjusted") +
  theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))
  
  
  # Extended Model ####
  
  # Importing results 
  res <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
  res[q.value < 0.001, q.value := round(q.value,3)]
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   ODI = res[exposure=="odi",])
  
  # filter MGS significant at the FDR p-value<0.05
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})
  
  # Create VennDiagram    
  venn.m3 <- ggvenn(mgs.fdr,
                    fill_color = c("orange","cornflowerblue","grey83"),
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .6, 
                    set_name_size = 4,
                    text_size = 3) +
    ggtitle("Extended adjustment") +
    theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))

  
  # Merge Venn diagrams
  venn <- plot_grid(venn.m1, venn.m2, venn.m3, labels = c("a","b","c"), label_size = 12, nrow=1,
                    label_y = .7)

  ggsave("Fig2.pdf",plot = venn,device = "pdf", path=output.plot)
  ggsave("Fig2.pdf",plot = venn,device = "pdf", path='/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
  
  

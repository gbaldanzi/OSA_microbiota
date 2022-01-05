# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-12-09

# Last update

# Venn diagram

# This code will produce a Venn Diagram from the correlations of AHI, BMI, 
# and T90 with MGSs in models 1 and 2

# Model 1 #### 

rm(list = ls())
# Loading packages 
pacman::p_load(data.table,ggplot2,ggvenn, tidyr, cowplot)

# input and output folders 
input = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing results 
res <- fread(paste0(input,"cor_all.var_mgs.tsv"))

res.list <- list(AHI = res[exposure=="ahi",],
                 T90 = res[exposure=="t90",],
                 ODI = res[exposure=="odi",],
                 BMI = res[exposure=="BMI",])


# filter MGS significant at the FDR p-value<0.05
res.list <- lapply(res.list, function(x){x[q.value>=0.001, q.value:=round(q.value, digits = 3)]})

mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})

# Create VennDiagram    
venn.m1 <- ggvenn(mgs.fdr,
                  fill_color = c("orange","cornflowerblue","grey83","green4"),
                  stroke_color = "white",
                  stroke_size = .2, 
                  show_percentage = F, 
                  fill_alpha = .6, 
                  set_name_size = 4,
                  text_size = 3) +
  ggtitle("Basic model") +
  theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))


# Model 2 ####

# Importing results 
res <- fread(paste0(input,"cor2_all.var_mgs.tsv"))

res.list <- list(AHI = res[exposure=="ahi",],
                 T90 = res[exposure=="t90",],
                 ODI = res[exposure=="odi",],
                 BMI = res[exposure=="BMI",])

# filter MGS significant at the FDR p-value<0.05
res.list <- lapply(res.list, function(x){x[q.value>=0.001, q.value:=round(q.value, digits = 3)]})

mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})

# Create VennDiagram    
venn.m2 <- ggvenn(mgs.fdr,
                  fill_color = c("orange","cornflowerblue","grey83","green4"),
                  stroke_color = "white",
                  stroke_size = .2, 
                  show_percentage = F, 
                  fill_alpha = .6, 
                  set_name_size = 4,
                  text_size = 3) +
  ggtitle("Fully adjusted model") +
  theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))

  
  # Merge Venn diagrams
  venn <- plot_grid(venn.m1, venn.m2, labels = c("A","B"), label_size = 12, nrow=1,
                    label_y = .7)

  ggsave("Venn_sleepapnea.png",plot = venn,device = "png", path=output.plot)
  ggsave("Venn_sleepapnea.pdf",plot = venn,device = "pdf", path='/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
  
  

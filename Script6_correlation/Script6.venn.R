# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2022-02-23

# Last update

# Venn diagram

# This code will produce a Venn Diagram from the correlations of AHI, ODI, 
# and T90 with MGSs in models 1 and 2

# Model 1 #### 

rm(list = ls())
# Loading packages 
pacman::p_load(data.table,ggplot2,ggvenn, tidyr, cowplot)

# results.folder and output folders 
results.folder <- "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing results 
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


# Model with BMI ####

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

  
  # Merge Venn diagrams
  venn <- plot_grid(venn.m1, venn.m2, labels = c("A","B"), label_size = 12, nrow=1,
                    label_y = .7)

  ggsave("Venn_sleepapnea.png",plot = venn,device = "png", path=output.plot)
  ggsave("Venn_sleepapnea.pdf",plot = venn,device = "pdf", path='/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
  
  

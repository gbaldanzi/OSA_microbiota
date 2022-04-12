# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-14

# Venn diagram

# This code will produce a Venn Diagram from the correlations of AHI imputed, ODI, 
# and T90 with gut microbitoa species in the basic and full model 

# Basic model #### 

rm(list = ls())
# Loading packages 
pacman::p_load(data.table,ggplot2,ggvenn, tidyr, cowplot)

# Function to adjust species names 
cutlast <- function(char,n){
  l <- nchar(char)
  a <- l-n+1
  return(substr(char,a,l))
}

# results.folder and output folders 
results.folder <- "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing results 
  res.t90.odi <- fread(paste0(results.folder,"cor_all.var_mgs.tsv"))
  res.t90.odi <- res.t90.odi[exposure %in% c("t90","odi"),]
  res.t90.odi[,MGS:=cutlast(MGS,9)]

  res.ahi <- fread(paste0(results.folder,"cor_ahi_imput_mgs.tsv"))
  res.ahi[,MGS:= gsub("_",".",MGS)]
  names(res.ahi) <- gsub("_",".",names(res.ahi))
  
  res <- rbind(res.ahi, res.t90.odi[, names(res.ahi),with=F])

  res[q.value >= 0.001, q.value := round(q.value,3)]
  
  # filter MGS significant at the FDR p-value<0.05
  res <- res[q.value < .05,]

  mgs.fdr <- list(AHI = res[exposure=="ahi",MGS],
                 T90 = res[exposure=="t90",MGS],
                 ODI = res[exposure=="odi",MGS])


# Create VennDiagram    
venn.m1 <- ggvenn(mgs.fdr,
                  fill_color = c("orange","cornflowerblue","grey83"),
                  stroke_color = "white",
                  stroke_size = .2, 
                  show_percentage = F, 
                  fill_alpha = .6, 
                  set_name_size = 4,
                  text_size = 3) +
  ggtitle("Basic model (AHI imputed)") +
  theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))


# Full model ####

# Importing results 
res.t90.odi <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
res.t90.odi <- res.t90.odi[exposure %in% c("t90","odi"),]
res.t90.odi[,MGS:=cutlast(MGS,9)]

res.ahi <- fread(paste0(results.folder,"cor2_ahi_imput_mgs.tsv"))
res.ahi[,MGS:= gsub("_",".",MGS)]
names(res.ahi) <- gsub("_",".",names(res.ahi))

res <- rbind(res.ahi, res.t90.odi[, names(res.ahi),with=F])

res[q.value >= 0.001, q.value := round(q.value,3)]

# filter MGS significant at the FDR p-value<0.05
res <- res[q.value < .05,]

mgs.fdr <- list(AHI = res[exposure=="ahi",MGS],
                T90 = res[exposure=="t90",MGS],
                ODI = res[exposure=="odi",MGS])

# Create VennDiagram    
venn.m2 <- ggvenn(mgs.fdr,
                  fill_color = c("orange","cornflowerblue","grey83"),
                  stroke_color = "white",
                  stroke_size = .2, 
                  show_percentage = F, 
                  fill_alpha = .6, 
                  set_name_size = 4,
                  text_size = 3) +
  ggtitle("Full model (AHI imputed)") +
  theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))

  
  # Merge Venn diagrams
  venn <- plot_grid(venn.m1, venn.m2, labels = c("A","B"), label_size = 12, nrow=1,
                    label_y = .7)

  ggsave("Venn_sleepapnea_imp.png",plot = venn,device = "png", path=output.plot)
  ggsave("Venn_sleepapnea_imp.pdf",plot = venn,device = "pdf", path='/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512')
  
  

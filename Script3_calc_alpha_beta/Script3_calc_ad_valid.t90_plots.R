# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will produce plots of the aitchison distance (euclidean distance of clr 
# transformed counts) between the participants with valid T90 measurements 


rm(list=ls())

# Loading packages
library(ape)
library(pheatmap)
library(data.table)
library(ggplot2)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")


# Import data
  valid.t90 <- readRDS("valid.t90_MGS.shannon_Upp.rds")

# Import Aitchision Distance matrix 
  AD <- fread(file = "/home/baldanzi/Datasets/sleep_SCAPIS/T90.ADmatrix.csv",sep=',')
  AD = as.matrix(AD)
  rownames(AD) = colnames(AD) 
  
# Load the PCoA results
  load('pc_AD_t90')
  
# Have datasets in the same row order 
  a <-  c("SCAPISid","t90cat")
  t90cat <- valid.t90[,a,with=F]
  t90cat <- t90cat[match(rownames(pcoa.ad$vectors),SCAPISid),]
  
  
# Plotting the first 2 principal components 
  
  # First: extracting the principal component value for every sample 
  dat.plot=data.frame(pcoa.ad$vectors, t90cat = t90cat$t90cat, SCAPISid=t90cat$SCAPISid)
  dat.plot$t90cat = factor(dat.plot$t90cat, levels = c("<10","10-20",">20"))
  
  # Second: creating the scatter plot 
  p1=ggplot(dat.plot,aes(x=Axis.1,y=Axis.2, color=t90cat))+
    geom_point(size=1.1) +
    stat_ellipse(type = "t", size=1.3) +
    ggtitle("Aitchison distance - T90%") +
    xlab(paste0("PCo1 \n (",round(100*pcoa.ad$values$Relative_eig[1],1),"% )")) +
    ylab(paste0("PCo2 \n (",round(100*pcoa.ad$values$Relative_eig[2],1),"% )")) +
    theme(axis.text.x = element_text(angle = 0),
          plot.title = element_text(hjust = 0.5))

ggsave("Atdist.T90.png", plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# T90 - Ait Dist Heatmap

paletteLength <- 500
myColor <- colorRampPalette(c("white","indianred3"))(paletteLength)
annotation = data.frame(t90cat = valid.t90$t90cat)
rownames(annotation) = colnames(AD)
annotation$t90cat = factor(annotation$t90cat, levels = c("<10","10-20",">20"))
rownames(AD) = colnames(AD)

# Heatmap
p1= pheatmap(AD,main="Heatmap",
             legend_labels = "Aitchison Distance - T90",
             color = myColor,
             show_rownames = FALSE, show_colnames = FALSE,
             cluster_rows = T,cluster_cols = T,
             annotation_row = annotation,
             annotation_col = annotation)

ggsave("T90.AD.heatmap.png", plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")
#----------------------------------------------------------------------------#

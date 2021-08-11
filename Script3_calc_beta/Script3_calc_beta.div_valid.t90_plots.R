# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will create plots and table describing beta
#diversity in relation to T90  in participants with valid T90

rm(list=ls())

# Loading packages
library(RColorBrewer)

# Import data
valid.t90 <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/valid.t90_MGS.shannon_Upp.rds")

# Import BC
BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/T90.BCmatrix.csv',sep=",")
BC = as.matrix(BC)
rownames(BC) = colnames(BC) 

# Import PCoA of BC 
load('/home/baldanzi/Datasets/sleep_SCAPIS/pc_BC_t90')

# Plotting the first 2 principal components 
# First: extracting the principal component value for every sample 
dat.plot=data.frame(pcoa.bray.t90$vectors, t90cat = valid.t90$t90cat)

# Second: creating the scatter plot ####
p4=ggplot(dat.plot, aes(x=Axis.1,y=Axis.2, color=t90cat))+
  geom_point(size=1.1) +
  stat_ellipse(type = "t", size=1.3) +
  ggtitle("Bray-curtis dissimilarity - T90") +
  xlab(paste0("PCo1 \n (",round(100*pcoa.bray.t90$values$Relative_eig[1],1),"% )")) +
  ylab(paste0("PCo2 \n (",round(100*pcoa.bray.t90$values$Relative_eig[2],1),"% )")) +
  theme(axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))


ggsave("BC.T90.png", plot = p4, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")


# BC heatmap - T90####

paletteLength <- 500
myColor <- colorRampPalette(c("white","indianred3"))(paletteLength)
annotation = data.frame(t90cat = valid.t90$t90cat)
rownames(annotation) = colnames(BC)
annotation$t90cat = factor(annotation$t90cat, levels = c("<10","10-20",">20"))

# Heatmap
p5= pheatmap(BC,main="Heatmap",
             legend_labels = "Bray-Curtis Dissimilarity - T90",
             color = myColor,
             show_rownames = FALSE, show_colnames = FALSE,
             cluster_rows = T,cluster_cols = T,
             annotation_row = annotation,
             annotation_col = annotation)

ggsave("T90.BC.heatmap.png", plot = p5, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

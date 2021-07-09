# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will create plots and table describing beta
#diversity in relation to AHI and OSA in participants with valid AHI

rm(list=ls())

# Loading packages
library(Hmisc)
library(compareGroups)
library(flextable)
library(pheatmap)
library(RColorBrewer)
library(grid)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")

# Import pheno data
valid.ahi <- readRDS("validsleep_MGS.shannon_Upp.rds")

# Import BC
BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv',sep=",")
BC = as.matrix(BC)
rownames(BC) = colnames(BC) 

# Import PCoA of BC 
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
load('pc_BC')


# Plotting the first 2 principal components 
    # First: extracting the principal component value for every sample 
      dat.plot=data.frame(pcoa.bray$vectors, OSAcat = valid.ahi$OSAcat)
      
    # Second: creating the scatter plot ####
      p4=ggplot(dat.plot, aes(x=Axis.1,y=Axis.2, color=OSAcat))+
          geom_point(size=1.1) +
          stat_ellipse(type = "t", size=1.3) +
          ggtitle("Bray-curtis dissimilarity - AHI") +
          xlab(paste0("PCo1 \n (",round(100*pcoa.bray$values$Relative_eig[1],1),"% )")) +
          ylab(paste0("PCo2 \n (",round(100*pcoa.bray$values$Relative_eig[2],1),"% )")) +
          theme(axis.text.x = element_text(angle = 0),
                plot.title = element_text(hjust = 0.5), face="bold")

  ggsave("BC.OSAcat.png", plot = p4, device = "png", 
         path = "/home/baldanzi/Sleep_apnea/Descriptive/")
  
  # Scatter plot - NoOSA vs Severe OSA ####
  severe_noosa <- valid.ahi$OSAcat[valid.ahi$OSA %in% c("no OSA","Severe")]

  dat.plot2=dat.plot[dat.plot$OSAcat %in% c("no OSA", "Severe"),]
  
  p4=ggplot(dat.plot2, aes(x=Axis.1,y=Axis.2, color=OSAcat))+
    geom_point(size=1.1) +
    stat_ellipse(type = "t", size=1.3) +
    ggtitle("Bray-curtis dissimilarity - AHI") +
    xlab(paste0("PCo1 \n (",round(100*pcoa.bray$values$Relative_eig[1],1),"% )")) +
    ylab(paste0("PCo2 \n (",round(100*pcoa.bray$values$Relative_eig[2],1),"% )")) +
    theme(axis.text.x = element_text(angle = 0),
          plot.title = element_text(hjust = 0.5), face="bold")
  
  ggsave("BC.OSAcat_noosa_severe.png", plot = p4, device = "png", 
         path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# BC heatmap - AHI ####

  paletteLength <- 500
  myColor <- colorRampPalette(c("white","indianred3"))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(0, 0.5, length.out=ceiling(paletteLength/2) + 1), 
                seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))

  annotation = data.frame(OSAcat = valid.ahi$OSAcat)
  rownames(annotation) = colnames(BC)

  # Heatmap
  p5= pheatmap(BC,main="Heatmap",
               legend_labels = "Bray-Curtis Dissimilarity",
               color = myColor,
               show_rownames = FALSE, show_colnames = FALSE,
               # breaks = myBreaks,
               cluster_rows = T,cluster_cols = T,
               annotation_row = annotation,
               annotation_col = annotation)

  ggsave("BC.heatmap.png", plot = p5, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")


# Box plots - BC ####
  BC = as.data.frame(BC)
  BC$SCAPISid = colnames(BC)

  # Merging with information on OSA cat 
  a = c("SCAPISid","OSAcat")
  BC2 = merge(BC, valid.ahi[,a,with=F], by="SCAPISid", all=T)
  BC2$OSAcat=factor(BC2$OSAcat, levels=c("no OSA","Mild","Moderate","Severe"))

  setDT(BC2)

  # Box plot - BC between "no OSA" and the other OSA groups
  a = c(BC2[OSAcat=="no OSA",SCAPISid],"SCAPISid","OSAcat")
  BCnoOSA= BC2[,a,with=F]
  BCnoOSA = BCnoOSA %>% gather(ID, distance, grep("5-", names(BCnoOSA), value = T))

  bp1 = ggplot(BCnoOSA, aes(x=OSAcat,y=distance,fill=OSAcat)) +
            geom_boxplot()+ xlab("Sleep apnea severity") + 
            ggtitle("BC between 'no OSA' and other categories")+
            theme(plot.title = element_text(hjust = .5 , size = 14,),
                  axis.title = element_text(size=14),
                  axis.text.x =  element_text(size=14),
                  legend.position = "none")

  ggsave('Boxplot_BC_noOSA.png', plot = bp1, path = "/home/baldanzi/Sleep_apnea/Descriptive/",
       device = "png")

# Box plot - BC  between "Severe" and the other OSA groups
  a = c(BC2[OSAcat=="Severe",SCAPISid],"SCAPISid","OSAcat")
  BCsevere= BC2[,a,with=F]
  BCsevere = BCsevere %>% gather(ID, distance, grep("5-", names(BCsevere), value = T))

  bp2 = ggplot(BCsevere, aes(x=OSAcat,y=distance,fill=OSAcat)) +
              geom_boxplot() + xlab("Sleep apnea severity") +
              ggtitle("BC between 'Severe' and other categories")+
              theme(plot.title = element_text(hjust = .5 , size = 14,),
                    axis.title = element_text(size=14),
                    axis.text.x =  element_text(size=14),
                    legend.position = "none")

  ggsave('Boxplot_BC_severe.png', plot = bp2, 
         path = "/home/baldanzi/Sleep_apnea/Descriptive/",
         device = "png")


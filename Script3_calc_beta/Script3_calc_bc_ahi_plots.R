# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will create plots and table describing beta
#diversity in relation to AHI and OSA in participants with valid AHI

  # Import BC
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv',sep=",")
  BC = as.matrix(BC)
  rownames(BC) = colnames(BC) 

  # Import PCoA of BC 
  load("/home/baldanzi/Datasets/sleep_SCAPIS/pc_BC")
  
  # pheno data 
  valid.ahi <- pheno[valid.ahi=='yes',]
  valid.ahi[OSAcat=="no OSA", OSAcat:="No Sleep Apnea"]
  valid.ahi[,OSAcat:=factor(OSAcat,levels = c("No Sleep Apnea","Mild","Moderate","Severe"))]

  # Plotting the first 2 principal components 
    # First: extracting the principal component value for every sample 
      dat.plot=data.frame(pcoa.bray$vectors, OSAcat = valid.ahi$OSAcat)
      
    # Second: creating the scatter plot ####
      
      centroids <- aggregate(cbind(Axis.1,Axis.2)~OSAcat,dat.plot,mean)
      f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
      se        <- aggregate(cbind(se.Axis.1=Axis.1,se.Axis.2=Axis.2)~OSAcat,dat.plot,f)
      centroids <- merge(centroids,se, by="OSAcat") 
      
      
      p4=ggplot(dat.plot, aes(x=Axis.1,y=Axis.2, color=OSAcat))+
         # geom_point(size=1.1) +
          geom_point(data=centroids, size=5)+
          geom_errorbar(data=centroids,aes(ymin=Axis.2-se.Axis.2,ymax=Axis.2+se.Axis.2),width=0.001)+
          geom_errorbarh(data=centroids,aes(xmin=Axis.1-se.Axis.1,xmax=Axis.1+se.Axis.1),height=0.0003) +
          #stat_ellipse(type = "t", size=1.3) +
         scale_color_manual(values=c("gray84","lightskyblue1","lightskyblue3","blue3")) +
          ggtitle("Bray-Curtis dissimilarity by sleep apnea severity") +
          xlab(paste0("PCo1 \n (",round(100*pcoa.bray$values$Relative_eig[1],1),"% )")) +
          ylab(paste0("PCo2 \n (",round(100*pcoa.bray$values$Relative_eig[2],1),"% )")) +
          xlim(-0.12, 0.04) + ylim(-0.08, .08) +
          theme_light()+
          theme(axis.text.x = element_text(angle = 0),
              plot.title = element_text(hjust = 0.5, face="bold", size=14),
              legend.text = element_text(size=10)) 
      
      p4$labels$colour <- "Sleep apnea severity"

  ggsave("BC.OSAcat.png", plot = p4, device = "png", 
         path = "/home/baldanzi/Sleep_apnea/Results/")
  
  ggsave("BC.OSAcat.pdf", plot = p4, device = "pdf", 
         path = "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/")
  
  
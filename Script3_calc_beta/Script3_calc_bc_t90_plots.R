# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# Last update: 2021-12-09

# This script will create plots and table describing beta
#diversity in relation to T90  in participants with valid T90

  # Import data
  valid.t90 <- pheno[valid.t90=='yes',]

  # Import BC
  BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/T90.BCmatrix.csv',sep=",")
  BC = as.matrix(BC)
  rownames(BC) = colnames(BC) 

  # Import PCoA of BC 
  load('/home/baldanzi/Datasets/sleep_SCAPIS/pc_BC_t90')

  # Plotting the first 2 principal components 
  # First: extracting the principal component value for every sample 
  dat.plot=data.frame(pcoa.bray.t90$vectors, t90cat = valid.t90$t90cat)
  
  dat.plot$t90cat <- factor(dat.plot$t90cat, 
                            levels = levels(dat.plot$t90cat) , 
                            labels = c("t1", "t2", "t3", "t4"))
  
  
  # Second: creating the scatter plot ####
  
  centroids <- aggregate(cbind(Axis.1,Axis.2)~t90cat,dat.plot,mean)
  f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
  se        <- aggregate(cbind(se.Axis.1=Axis.1,se.Axis.2=Axis.2)~t90cat,dat.plot,f)
  centroids <- merge(centroids,se, by="t90cat") 
  
  
  p2=ggplot(dat.plot, aes(x=Axis.1,y=Axis.2, color=t90cat))+
    # geom_point(size=1.1) +
    geom_point(data=centroids, size=4)+
    geom_errorbar(data=centroids,aes(ymin=Axis.2-se.Axis.2,ymax=Axis.2+se.Axis.2),width=0.001)+
    geom_errorbarh(data=centroids,aes(xmin=Axis.1-se.Axis.1,xmax=Axis.1+se.Axis.1),height=0.0003) +
    #stat_ellipse(type = "t", size=1.3) +
    scale_color_manual(values=c("gray84","lightskyblue1","lightskyblue3","blue3")) +
    ggtitle("T90 severity groups") +
    xlab(paste0("PCo1 \n (",round(100*pcoa.bray.t90$values$Relative_eig[1],1),"% )")) +
    ylab(paste0("PCo2 \n (",round(100*pcoa.bray.t90$values$Relative_eig[2],1),"% )")) +
    xlim(-0.08, 0.08) + ylim(-0.08, .08) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 0),
          axis.title = element_text(size=7) , 
          axis.text = element_text(size=6) , 
          plot.title = element_text(hjust = 0.5, face="bold", size=10),
          legend.text = element_text(size=8),
          legend.position = "bottom") 
  
  p2$labels$colour <- "T90 severity groups"
  
  p2 <- p2 + theme(legend.title = element_blank()) 
  
  p2 <- p2 + guides(color=guide_legend(nrow=2))
  
  ggsave("BC.t90cat.png", plot = p2, device = "png", 
         path = "/home/baldanzi/Sleep_apnea/Results/")
  
  ggsave("BC.t90cat.pdf", plot = p2, device = "pdf", 
         path = "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/")
  
  





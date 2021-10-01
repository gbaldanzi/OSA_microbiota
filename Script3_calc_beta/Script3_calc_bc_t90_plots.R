# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

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

  # Second: creating the scatter plot ####
  p4=ggplot(dat.plot, aes(x=Axis.1,y=Axis.2, color=t90cat))+
      geom_point(size=1.1) +
      stat_ellipse(type = "t", size=1.3) +
      ggtitle("Bray-curtis dissimilarity - T90") +
      xlab(paste0("PCo1 \n (",round(100*pcoa.bray.t90$values$Relative_eig[1],1),"% )")) +
      ylab(paste0("PCo2 \n (",round(100*pcoa.bray.t90$values$Relative_eig[2],1),"% )")) +
      theme(axis.text.x = element_text(angle = 0),
            plot.title = element_text(hjust = 0.5, face='bold'))


  ggsave("BC.T90.png", plot = p4, device = "png", 
        path = "/home/baldanzi/Sleep_apnea/Descriptive/")



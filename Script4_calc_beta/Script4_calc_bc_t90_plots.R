# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# # PCoA separated by T90 severity groups 

  # Import data
  valid.t90 <- pheno[valid.t90=='yes',]
  
  valid.t90[,t90cat:= factor(t90cat, 
                            levels = levels(t90cat) , 
                            labels = c("T90=0", "t1", "t2", "t3")) ]

  # Import BC
  BC <-  fread(paste0(input,'T90.BCmatrix.csv'),sep=",")
  BC <-  as.matrix(BC)
  rownames(BC) <-  colnames(BC) 

  # Import PCoA of BC 
  load(paste0(input,'pc_BC_t90'))

  # Plotting the first 2 principal coordinate 
  # First: extracting the principal coordinate values for every sample 
  dat.plot=data.frame(pcoa.bray.t90$vectors, t90cat = valid.t90$t90cat)
  
  dat.plot$t90cat <- factor(dat.plot$t90cat, 
                            levels = levels(dat.plot$t90cat) , 
                            labels = c("T90 = 0", "t1", "t2", "t3"))
  
  
  # Second: creating the scatter plot ####
  
  centroids <- aggregate(cbind(Axis.1,Axis.2)~t90cat,dat.plot,mean)
  f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
  se        <- aggregate(cbind(se.Axis.1=Axis.1,se.Axis.2=Axis.2)~t90cat,dat.plot,f)
  centroids <- merge(centroids,se, by="t90cat") 
  
  Ns <- with(valid.t90, table(t90cat))
  
  p2=ggplot(dat.plot, aes(x=Axis.1,y=Axis.2, color=t90cat))+
    geom_point(data=centroids, size=2)+
    geom_errorbar(data=centroids,aes(ymin=Axis.2-se.Axis.2,ymax=Axis.2+se.Axis.2),width=0.001)+
    geom_errorbarh(data=centroids,aes(xmin=Axis.1-se.Axis.1,xmax=Axis.1+se.Axis.1),height=0.0003) +
    scale_color_manual(labels = paste0(names(Ns), " (n = ",formatC(Ns,format = "d", big.mark = ","),")"),
                       values=c("gray84","lightskyblue1","lightskyblue3","blue3")) +
    ggtitle("T90 severity groups") +
    xlab(paste0("PCo1 \n (",round(100*pcoa.bray.t90$values$Relative_eig[1],1),"%)")) +
    ylab(paste0("PCo2 \n (",round(100*pcoa.bray.t90$values$Relative_eig[2],1),"%)")) +
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
  
  p2 <- p2 + guides(color=guide_legend(nrow=4))
  

  





# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# PCoA separated by ODI severity groups 

  # PCoA ####
  
  centroids <- aggregate(cbind(Axis.1,Axis.2)~odicat,dat.plot,mean)
  f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
  se        <- aggregate(cbind(se.Axis.1=Axis.1,se.Axis.2=Axis.2)~odicat,dat.plot,f)
  centroids <- merge(centroids,se, by="odicat") 
  
  Ns <- with(valid.t90, table(odicat))
  
  p3=ggplot(dat.plot, aes(x=Axis.1,y=Axis.2, color=odicat))+
    geom_point(data=centroids, size=2)+
    geom_errorbar(data=centroids,aes(ymin=Axis.2-se.Axis.2,ymax=Axis.2+se.Axis.2),width=0.0003)+
    geom_errorbarh(data=centroids,aes(xmin=Axis.1-se.Axis.1,xmax=Axis.1+se.Axis.1),height=0.0003) +
    scale_color_manual(labels = paste0(names(Ns), " (n = ",Ns,")"),
                       values=c("gray84","lightskyblue1","lightskyblue3","blue3")) +
    ggtitle("ODI severity groups") +
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
  
  p3$labels$colour <- "ODI severity groups"
  
  p3 <- p3 + theme(legend.title = element_blank()) 
  
  p3 <- p3 + guides(color=guide_legend(nrow=4))
  



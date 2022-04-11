# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

library(data.table)
library(tidyverse)
library(cowplot)

# PCoA separated by AHI/T90/ODI severity groups 

  input <- "/home/baldanzi/Datasets/sleep_SCAPIS/"
  # Import data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  pheno[,t90cat := factor(t90cat, levels = levels(t90cat), 
                          labels = c("T90 = 0", "t1", "t2","t3"))]

  # Import BC
  BC <-  fread(paste0(input,'BCmatrix.tsv'))
  BC_rownames <- BC$rownames
  BC[,rownames:=NULL]
  BC <-  as.matrix(BC)
  colnames(BC) <- rownames(BC) <- BC_rownames

  # Import PCoA of BC 
  pcoa.bray.ahi <- readRDS(paste0(input,"pcoa_BC_AHI.rds"))
  pcoa.bray.t90 <- readRDS(paste0(input,"pcoa_BC_T90.ODI.rds"))
  
  # Improve visualization by rotating 180°
  pcoa.bray.t90$vectors[,1] <- pcoa.bray.t90$vectors[,1]*(-1)
  
  # PCOA.plot function ####
  pcoa.plot <- function(data,groups,pcoa) {
    
  df <- data[SCAPISid %in% rownames(pcoa$vectors),]

  # Plotting the first 2 principal coordinates 
    # First: extracting the principal coordinate values for every sample 
      dat.plot=data.frame(pcoa$vectors, df[,groups,with=F])
      names(dat.plot) <- c(colnames(pcoa$vectors),"cat")
      
      
    # Second: creating the scatter plot ####
      
      centroids <- aggregate(cbind(Axis.1,Axis.2)~cat,dat.plot,mean)
      f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
      se        <- aggregate(cbind(se.Axis.1=Axis.1,se.Axis.2=Axis.2)~cat,dat.plot,f)
      centroids <- merge(centroids,se, by="cat") 
      
      Ns <- table(dat.plot$cat)
      
      .limits <- range(centroids$Axis.1)
      .limits[1] <- .limits[1]-.01
      .limits[2] <- .limits[2]+.02
      
      p1=ggplot(centroids, aes(x=Axis.1,y=Axis.2, color=cat)) +
          geom_point(size=2) +
          geom_errorbar(aes(ymin=Axis.2-se.Axis.2,ymax=Axis.2+se.Axis.2),width=0.001) +
          geom_errorbarh(aes(xmin=Axis.1-se.Axis.1,xmax=Axis.1+se.Axis.1),height=0.0003) +
         scale_color_manual(labels =paste0(names(Ns), " (n = ",formatC(Ns,format = "d",big.mark = ","),")") , 
           values=c("gray84","lightskyblue1","lightskyblue3","blue3")) +
          ggtitle("XXX severity groups") +
          xlab(paste0("PCo1 \n (",round(100*pcoa$values$Relative_eig[1],1),"%)")) +
          ylab(paste0("PCo2 \n (",round(100*pcoa$values$Relative_eig[2],1),"%)")) +
          scale_x_continuous(breaks = seq(-0.06,0.09,by=.03), limits = .limits ) + 
          scale_y_continuous(breaks = seq(-0.06,0.09,by=.03), limits = .limits ) +
          theme_light()+
          theme(axis.text.x = element_text(angle = 0),
                axis.title = element_text(size=7) ,
                axis.text = element_text(size=6) , 
              plot.title = element_text(hjust = 0.5, face="bold", size=10),
              legend.text = element_text(size=8),
              legend.position = "bottom") 
      
      p1$labels$colour <- "Sleep apnea severity"

      p1 <- p1 + theme(legend.title = element_blank())  
      
      p1 <- p1 + guides(color=guide_legend(nrow = 4))
      return(p1)
  }
  
  
  p.ahi <- pcoa.plot(pheno,"OSAcat",pcoa.bray.ahi)
  p.ahi$labels$title <- "AHI severity groups"
  
  p.t90 <- pcoa.plot(pheno,"t90cat",pcoa.bray.t90)
  p.t90$labels$title <- "T90 severity groups"
  
  p.odi <- pcoa.plot(pheno,"odicat",pcoa.bray.t90)
  p.odi$labels$title <- "ODI severity groups"
  
  pcoa.plot.merged <- plot_grid(NULL,NULL,NULL,p.ahi,p.t90,p.odi,NULL,NULL,NULL,labels = c("","","","a","b","c"), label_size = 12, nrow=3,
                                label_y = 0.97, ncol=3, rel_heights = c(.5,1,.4))
  ggsave("PCoA_sleepapnea.pdf", plot = pcoa.plot.merged, device = "pdf", 
         path = '/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/')
  
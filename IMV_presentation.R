# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

library(data.table)
library(tidyverse)
library(cowplot)
library(RColorBrewer)

# PCoA separated by AHI/T90/ODI severity groups 

input <- "/home/baldanzi/Datasets/sleep_SCAPIS/"
# Import data
pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))

# Import BC
BC <-  fread(paste0(input,'BCmatrix.tsv'))
BC_rownames <- BC$rownames
BC[,rownames:=NULL]
BC <-  as.matrix(BC)
colnames(BC) <- rownames(BC) <- BC_rownames

# Import PCoA of BC 
pcoa <- readRDS(paste0(input,"pcoa_BC_AHI.rds"))


df <- pheno[SCAPISid %in% rownames(pcoa$vectors),]

df[,OSAcat:=factor(OSAcat, levels = c("No OSA","Mild","Moderate","Severe"),
                   labels = c("AHI<5","AHI 5-14.9","AHI 15-29.9","AHI>30"))]

# Plotting the first 2 principal coordinates 
# First: extracting the principal coordinate values for every sample 
dat.plot=data.frame(pcoa$vectors, df[,OSAcat])
names(dat.plot) <- c(colnames(pcoa$vectors),"cat")


# Second: creating the scatter plot ####

centroids <- aggregate(cbind(Axis.1,Axis.2)~cat,dat.plot,mean)
f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
se        <- aggregate(cbind(se.Axis.1=Axis.1,se.Axis.2=Axis.2)~cat,dat.plot,f)
centroids <- merge(centroids,se, by="cat") 

Ns <- table(dat.plot$cat)

.limits <- c(-.4,0.48)


p1=ggplot(centroids, aes(x=Axis.1,y=Axis.2, color=cat)) +
  #geom_point(size=2) +
  geom_point(data=dat.plot, aes(x=Axis.1,y=Axis.2, color=cat),alpha=.7) +
  #geom_errorbar(aes(ymin=Axis.2-se.Axis.2,ymax=Axis.2+se.Axis.2),width=0.001) +
  #geom_errorbarh(aes(xmin=Axis.1-se.Axis.1,xmax=Axis.1+se.Axis.1),height=0.0003) +
  scale_color_manual(labels = names(Ns), 
                     values=c("#FFFFCC","#A1DAB4","#41B6C4","blue4")) +
  xlab(paste0("PCo1 \n (",round(100*pcoa$values$Relative_eig[1],1),"%)")) +
  ylab(paste0("PCo2 \n (",round(100*pcoa$values$Relative_eig[2],1),"%)")) +
  scale_x_continuous(breaks = seq(-0.25,0.5,by=0.25), limits = .limits ) + 
  scale_y_continuous(breaks = seq(-0.25,0.5,by=0.25), limits = .limits ) +
  #theme_light()+
  theme(axis.text.x = element_text(angle = 0),
        axis.title = element_text(size=10) ,
        axis.text = element_text(size=8) , 
        plot.title = element_text(hjust = 0.5, face="bold", size=10),
        legend.text = element_text(size=12),
        legend.position = "right") 

p1$labels$colour <- "Sleep apnea severity"

p1 <- p1 + theme(legend.title = element_blank())  

p1 <- p1 + guides(color=guide_legend(nrow = 4))

ggsave("legend.fig1.png",p1)

p1 <- p1 + theme(legend.position = "none")


ggsave("presentation.fig1.png",p1)

# Plot 2 ####

p2=ggplot(centroids, aes(x=Axis.1,y=Axis.2, color=cat)) +
  geom_point(size=3.3, shape = 15) +
  geom_point(data=dat.plot, aes(x=Axis.1,y=Axis.2, color=cat),alpha=.15) +
  geom_errorbar(aes(ymin=Axis.2-se.Axis.2,ymax=Axis.2+se.Axis.2),width=0.00) +
  geom_errorbarh(aes(xmin=Axis.1-se.Axis.1,xmax=Axis.1+se.Axis.1),height=0.000) +
  scale_color_manual(values=c("#FFFFCC","#A1DAB4","#41B6C4","blue4")) +
 
  stat_ellipse(data=dat.plot, aes(x=Axis.1,y=Axis.2, color=cat),type = "t", 
                size=1.0) +
  xlab(paste0("PCo1 \n (",round(100*pcoa$values$Relative_eig[1],1),"%)")) +
  ylab(paste0("PCo2 \n (",round(100*pcoa$values$Relative_eig[2],1),"%)")) +
  scale_x_continuous(breaks = seq(-0.25,0.5,by=0.25), limits = .limits ) + 
  scale_y_continuous(breaks = seq(-0.25,0.5,by=0.25), limits = .limits ) +
  #theme_light()+
  theme(axis.text.x = element_text(angle = 0),
        axis.title = element_text(size=10) ,
        axis.text = element_text(size=8) , 
        plot.title = element_text(hjust = 0.5, face="bold", size=10),
        #legend.text = element_text(size=12),
        legend.position = "none") 

p2$labels$colour <- "Sleep apnea severity"

p2 <- p2 + theme(legend.title = element_blank())  

p2 <- p2 + guides(color=guide_legend(nrow = 4))


ggsave("presentation.fig2.png",p2)


# Enrichment analysis 
  results.folder <-  "/home/baldanzi/Sleep_apnea/Results/"
  res.gmm = fread(paste0(results.folder,"ea_GMM_pos_new.tsv"))
  res.gmm <- res.gmm[exposure=="t90",]
  res.gmm[,sig := ifelse(q.value<.05,"sig","non.sig")]
  res.gmm[,id:= rank(pval)]
  
  gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  setnames(gmm.names,"Module","modules")
  gmm.names[,Name:=str_to_title(Name)]
  gmm.names[,Name:=gsub("Ii","II",Name)]
  
  res.gmm <- merge(res.gmm,gmm.names,by.x="pathway",by.y="modules",all.x=T,all.y=F)
  
  
  p3 <- ggplot(res.gmm,aes(y=NES,x=-log(pval),col=sig)) + 
    geom_point() + 
    scale_color_manual(values=c("gray84","blue3")) + 
    xlab("-log(p-value)") + ylab("Normalized enrichment score") + 
    geom_vline(xintercept=-log(0.006), linetype = "dashed", size=1,color="darkred") +
    ggtitle("Enriched pathways") +
    theme_light() + 
    geom_text(aes(label=ifelse(q.value<0.05,as.character(id),'')),hjust=0,vjust=0)

  p3$theme$plot.title <- element_text(size=18,hjust=0.5)
  
  p3$theme$axis.text <- element_blank()
  
  p3$theme$axis.title = element_text(size=12)

  p3$theme$legend.position = "none"
   
  ggsave("figure.ea.png",plot=p3)
    






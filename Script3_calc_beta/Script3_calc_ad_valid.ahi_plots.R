# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will produce plots of the aitchison distance (euclidean distance of clr 
# transformed counts) between the participants with valid ahi 


rm(list=ls())
#Import data 

# Import pheno data
  valid.ahi <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/validsleep_MGS.shannon_Upp.rds")

# Import Aitchison distance matrix 
  AD = fread("/home/baldanzi/Datasets/sleep_SCAPIS/OSA.aitchison_distmatrix.csv")
  AD = as.matrix(AD)
  rownames(AD) = colnames(AD) 

# Load PCoA results 
  load("/home/baldanzi/Datasets/sleep_SCAPIS/pc_AD")
  
# Have datasets in the same row order 
  a <-  c("SCAPISid","OSAcat")
  OSAcat <- valid.ahi[,a,with=F]
  OSAcat <- OSAcat[match(row.names(pcoa.ad$vectors),SCAPISid),]

# Plotting the first 2 principal components 
# First: extracting the principal component value for every sample 
dat.plot=data.frame(pcoa.ad$vectors, OSAcat = OSAcat$OSAcat, SCAPISid=OSAcat$SCAPISid)
dat.plot$OSAcat = factor(dat.plot$OSAcat, levels=c("no OSA", "Mild", "Moderate", "Severe"))

# Second: creating the scatter plot 
p1=ggplot(dat.plot,aes(x=Axis.1,y=Axis.2, color=OSAcat))+
  geom_point(size=1.1)+
  stat_ellipse(type = "t", size=1.3) +
  ggtitle("Aitchison distance - OSA severity") +
  xlab(paste0("PCo1 \n (",round(100*pcoa.ad$values$Relative_eig[1],1),"% )")) +
  ylab(paste0("PCo2 \n (",round(100*pcoa.ad$values$Relative_eig[2],1),"% )")) +
  theme(axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))

ggsave("Atdist.OSAcat.png", plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# Creating a Aitchison distance heatmap ####

paletteLength <- 500
myColor <- colorRampPalette(c("white","indianred3"))(paletteLength)

#annotation 
annotation = data.frame(OSAcat = OSAcat$OSAcat)
rownames(annotation) = colnames(AD)
annotation$OSAcat = factor(annotation$OSAcat, levels=c("no OSA", "Mild", "Moderate", "Severe"))

# Heatmap
p2=pheatmap(AD,main="Heatmap",
            legend_labels = "Aitchison Distance",
            color = myColor,show_rownames = FALSE, show_colnames = FALSE,
            cluster_rows = T,cluster_cols = T, 
            annotation_row = annotation,
            annotation_col = annotation)

ggsave("AitDist.heatmap.png", plot = p2, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# Merge clr-data with valid sleep recordings 
#a = c("SCAPISid", "OSAcat", "ahi", "age", "Sex")
#count.clr=merge(count.clr,valid.ahi[,a,with=F], by="SCAPISid",all=F, all.x=F, all.y=T)

# Saving clr-transformed data 
#fwrite(count.clr,"/home/baldanzi/Datasets/sleep_SCAPIS/validsleep.mgs_clr_count.tsv",sep='\t')

#Saving the atchison distances 
#atdist_df = as.data.frame(atdist)
#atdist_df$SCAPISid = count$SCAPISid
#atdist_df$OSAcat = count$OSAcat
#fwrite(atdist_df,"/home/baldanzi/Datasets/sleep_SCAPIS/OSA.aitchison_distmatrix.csv",sep=',')

# Saving the principal components analysis 
#save(pc.atdist, file = '/home/baldanzi/Datasets/sleep_SCAPIS/pc_aitdist')


# Box plots - Aitchison distances ####

# Box plot - Ait. Dist. between no OSA and the other categories
  OSAcat <- OSAcat[match(row.names(AD),SCAPISid),]
  AD = data.frame(AD)
  AD <- cbind(AD, OSAcat$OSAcat, OSAcat$SCAPISid)
  names(AD) = c(rownames(AD),"OSAcat","SCAPISid") 
  a = c(AD[AD$OSAcat=="no OSA","SCAPISid"],"SCAPISid","OSAcat")
  ADnoOSA= AD[,a]
  ADnoOSA = ADnoOSA %>% gather(ID, distance, grep("5-", names(ADnoOSA), value = T))

bp1 = ggplot(ADnoOSA, aes(x=OSAcat,y=distance,fill=OSAcat)) +
  geom_boxplot()+ xlab("Sleep apnea severity") + 
  ggtitle("Ait. Dist. between 'no OSA' and other categories")+
  theme(plot.title = element_text(hjust = .5 , size = 14,),
        axis.title = element_text(size=14),
        axis.text.x =  element_text(size=14),
        legend.position = "none")

ggsave('Boxplot_AD_noOSA.png', plot = bp1, path = "/home/baldanzi/Sleep_apnea/Descriptive/",
       device = "png")

# Box plot - Ait. Dist. between Severe and the other categories
  a = c(AD[AD$OSAcat=="Severe","SCAPISid"],"SCAPISid","OSAcat")
  ADsevere= AD[,a]
  ADsevere = ADsevere %>% gather(ID, distance, grep("5-", names(ADsevere), value = T))

bp2 = ggplot(ADsevere, aes(x=OSAcat,y=distance,fill=OSAcat)) +
  geom_boxplot() + xlab("Sleep apnea severity") +
  ggtitle("Ait. Dist. between 'Severe' and other categories")+
  theme(plot.title = element_text(hjust = .5 , size = 14,),
        axis.title = element_text(size=14),
        axis.text.x =  element_text(size=14),
        legend.position = "none")

ggsave('Boxplot_AD_severe.png', plot = bp2, path = "/home/baldanzi/Sleep_apnea/Descriptive/",
       device = "png")


# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics continued : script to DESCRIBE alpha and beta-diversity

# This script will calculate alpha diversity and create plots and table describing alpha
#diversity in relation AHI and OSA, as well as other covarites. 

# Descriptive Statistics 

rm(list=ls())

# Loading packages
library(tidyverse)
library(data.table)
library(vegan)
library(Hmisc)
library(compareGroups)
library(flextable)
library(pheatmap)
library(RColorBrewer)
library(robCompositions)
library(grid)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")

# Uploading dataset with variables of class "factor"
load("data_table1")
# Importing data on individuals with valid.ahi measurement 
valid.ahi = fread("validsleep.MGS.Upp.tsv", header=T, na.strings=c("", "NA"))
#merging with the data set that has variables as factor 
a = names(dat1[,-"SCAPISid"])
valid.ahi[,(a):=NULL]
valid.ahi = merge(dat1,valid.ahi, by = "SCAPISid", all=T, all.x=F, all.y=F)


#Shannon diversity ####
a = grep("____",names(valid.ahi),value=T) # vector with MGS names 
valid.ahi$shannon=diversity(valid.ahi[,a, with=F],index="shannon") #estimating shannon per individual

#fwrite(valid.ahi, file = "validsleep_MGS.shannon.BC_Upp.tsv", sep = "\t")

summary(valid.ahi[,shannon])
valid.ahi[,mean(shannon), by=OSAcat]
valid.ahi[,median(shannon), by=OSAcat]

# Compare the Shannon index across groups using Kruskal-Wallis #### 
res=with(valid.ahi, kruskal.test(shannon ~ OSAcat))

# Scatter plot: Shannon index against AHI 
p1 = ggplot(data=valid.ahi) + 
  geom_point(aes(x=ahi, y= shannon)) + ggtitle("Shannon index and AHI") + 
  xlab("Apnea-hypopnea index") +
  ylab("Shannon index") +
  geom_vline(xintercept = 5, color = "red", linetype = "twodash") + 
  geom_vline(xintercept = 15, color = "red", linetype = "twodash") + 
  geom_vline(xintercept = 30, color = "red", linetype = "twodash")

ggsave("scatter.shannon.ahi", plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

p2 = ggplot(data=valid.ahi) + 
  geom_point(aes(x=rank(ahi), y= rank(shannon))) + ggtitle("Rank Shannon index and rank AHI") + 
  xlab("Apnea-hypopnea index") +
  ylab("Shannon index") 

ggsave("rank_scatter.shannon.ahi", plot = p2, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

p3=ggplot(valid.ahi, aes(x=OSAcat, y=shannon)) +
  geom_violin(trim=FALSE, fill="indianred3")+
  geom_boxplot(width=0.1)+
  ggtitle(paste0("Shannon index by OSA severity - Kruskal-Wallis: p.value=",formatC(res$p.value,format= "e", digits=2)))+
  scale_x_discrete(labels=paste(names(table(valid.ahi$OSAcat))," (n=",table(valid.ahi$OSAcat),")",sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("OSA severity")+
  ylab("Shannon Diversity Index")

ggsave("violin.shannon.OSAcat", plot = p3, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

p3_1=ggplot(data=valid.ahi,aes(x=shannon))+
  geom_histogram()+
  ggtitle("Histogram - Shannon diversity")
ggsave("hist.shannon", plot = p3_1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# Bray-curtis dissimilarity #### 

# Create dataset containing exclusively the MGS as columns/variables 
a=grep("____",names(valid.ahi),value=T)
MGS=valid.ahi[,a, with=F]

# Estimating the BC index - creates a matrix with the BC index between individual samples 
#BC=as.matrix(vegdist(MGS,method="bray"))
#rownames(BC)=colnames(BC)=valid.ahi$SCAPISid
BC = fread('/home/baldanzi/Datasets/sleep_SCAPIS/OSA.BCmatrix.csv',sep=",")
BC = as.matrix(BC)
rownames(BC) = colnames(BC) 

# Principal coordinates analysis based on BC 
#pc.bray=prcomp(BC) 
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
load('pc_BC')

# Extracting the proportion of variance explained for each principal component (eigenvalues)
prop=data.frame("PC"=names(summary(pc.bray)$importance[2,]),
                "value"=summary(pc.bray)$importance[2,]*100,
                "cumulative_prop" = summary(pc.bray)$importance[3,]*100,
                stringsAsFactors = T)

# Saving data with shannon index and principal components of BC 
# Saving datasets ####
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
#save(pc.bray, file = 'pc_BC')
#fwrite(valid.ahi, file = "validsleep_MGS.shannon.BC_Upp.tsv", sep = "\t")
#fwrite(BC, file = "OSA.BCmatrix.csv", sep = ",")

# Plotting the first 2 principal components 
# First: extracting the principal component value for every sample 
dat.plot=data.frame(pc.bray$x, OSAcat = valid.ahi$OSAcat)

# Second: creating the scatter plot ####
p4=ggplot(dat.plot,aes(x=PC1,y=PC2, color=OSAcat))+
  geom_point(size=1.1)+
  stat_ellipse(type = "t", size=1.3) +
  ggtitle("Bray-curtis dissimilarity") +
  xlab(paste0("PC1 \n (",round(summary(pc.bray)$importance[2,1]*100,2),"% )")) +
  ylab(paste0("PC2 \n (",round(summary(pc.bray)$importance[2,2]*100,2),"% )")) +
  theme(axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))

ggsave("BC.OSAcat.png", plot = p4, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# BC heatmap ####

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

# Box plot - BC disimilarity between no OSA and the other categories
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

# Box plot - BC disimilarity between Severe and the other categories
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

ggsave('Boxplot_BC_severe.png', plot = bp2, path = "/home/baldanzi/Sleep_apnea/Descriptive/",
       device = "png")

#-------------------------------------------------------------------------#
# MGS prevalence ####

# MGS prevalence is the number of participants in which a certain MGS is present

# presence-absence transformation:
# If a species is present, it is transformed to 1
# If it is absent, it is transformed to zero
a = grep("____H",names(valid.ahi))
data_pa <- decostand(x = valid.ahi[,a,with=F], "pa")

# calculate sum per species
data_sum <- apply(data_pa, 2, sum)

p1 = ggplot(data = as.data.frame(data_sum), aes(x=data_sum)) + 
  geom_histogram(bins = 70,color="black",fill="white")  +
  ggtitle("Histogram MGS prevalence") +  
  xlab(paste("Present in at least 100 indiv: ", 
             length(data_sum[data_sum>=100])," (",
             round(length(data_sum[data_sum>=100])/length(data_sum),2)*100,"%)\n",
             "Present in at least 50 indiv: ",
             length(data_sum[data_sum>=50])," (",
             round(length(data_sum[data_sum>=50])/length(data_sum),2)*100,"%)\n",
             sep="")) + 
  geom_vline(xintercept = 50, color = "black", linetype = "twodash") + 
  geom_vline(xintercept = 100, color = "black", linetype = "twodash") +
  theme(plot.title = element_text(hjust = .5,size = 14,face = "bold"))

ggsave("hist.MGS.prevalence",plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

#--------------------------------------------------------------------------#
# AITCHISON DISTANCE ####
# Beta diversity based on Aitchison distance. ####
print("Aitchison distance")


# Aitchison distance the Eucledian distance between CLR transformed count data. 

# Importing count data:
#count = readRDS('/home/baldanzi/Datasets/MGS/upugut03.mgsCounts.rds')

# Summary of total counts per sample 
#row_sum=rowSums(count)
#summary(row_sum)
#  Min. 1st Qu.  Median    Mean  3rd Qu.    Max. 
# 1397  587162  657876   672899  740013   2016207 

# Number of zeros in the data
#sum(count==0)
# 8142757

# Zeros are replaced with 1 (there were no previous 1 in the data)
#count[count==0]=1

# Central log-transformation of count data 
#count.clr = cenLR(count)
#count=as.data.frame(count)
#count.clr = as.data.frame(count.clr$x.clr)


#To fix some added characters to the names
#rownames(count)[grep("gut.+_b",rownames(count),value=F)]=gsub("_b","",grep("gut.+_b",rownames(count),value=T))
#rownames(count)[grep("gut.+b",rownames(count),value=F)]=gsub("b","",grep("gut.+b",rownames(count),value=T))
#rownames(count)[grep("gut.+c",rownames(count),value=F)]=gsub("c","",grep("gut.+c",rownames(count),value=T))

#rownames(count.clr)[grep("gut.+_b",rownames(count.clr),value=F)]=gsub("_b","",grep("gut.+_b",rownames(count.clr),value=T))
#rownames(count.clr)[grep("gut.+b",rownames(count.clr),value=F)]=gsub("b","",grep("gut.+b",rownames(count.clr),value=T))
#rownames(count.clr)[grep("gut.+c",rownames(count.clr),value=F)]=gsub("c","",grep("gut.+c",rownames(count.clr),value=T))


# Adjusting the ID variable 
#count$sample.id=rownames(count)
#count.clr$sample.id=rownames(count.clr)



#Merging with the SCAPISid
#selected_id = fread("/home/baldanzi/Datasets/scapis_idkey/selected_participants_id_sample_id.tsv", 
#                  header=T, na.strings=c("", "NA"))
#count = merge(count, selected_id, by="sample.id",all=F,all.x=T, all.y=F)
#count.clr = merge(count.clr, selected_id, by="sample.id",all=T,all.x=F, all.y=F)

#ids=fread("/home/baldanzi/Datasets/scapis_idkey/id_conversion.txt", header=T) 
#ids$id = as.character(ids$export_id)
#setnames(ids, "subject_id", "SCAPISid")
#count = merge(count, ids, by = "id", all=F, all.x=T,all.y=F)
#count.clr = merge(count.clr, ids, by = "id", all=F, all.x=T,all.y=F)

# Only keeping the individuals that have valid sleep recording
#count=count[count$SCAPISid %in% valid.ahi[,SCAPISid],]

# Merging counts with valid.ahi 
#a = c("SCAPISid", "OSAcat", "ahi", "age", "Sex")
#count=merge(count,valid.ahi[,a,with=F], by="SCAPISid",all=F, all.x=F, all.y=T)

# Aitchison distance ####
#rownames(count)=count$SCAPISid
#mgs_count = as.matrix(count[,grep("HG3A",names(count))])


#t0 = Sys.time()
#atdist=aDist(mgs_count)
#t1 = Sys.time()
#print(paste0("time=" ,t1-t0))

# Transform aitchison distances into a matrix
#atdist_matrix=as.matrix(atdist)

# REMOVE THIS ####
atdist_df = fread("/home/baldanzi/Datasets/sleep_SCAPIS/OSA.aitchison_distmatrix.csv",sep=',')
rownames(atdist_df)=atdist_df$SCAPISid
a = grep("-",names(atdist_df), value=T)
atdist_matrix = as.matrix(atdist_df[,a,with=F])

rownames(atdist_matrix) = colnames(atdist_matrix)

count = fread("/home/baldanzi/Datasets/sleep_SCAPIS/validsleep.mgs_clr_count.tsv",sep='\t')
#####

# PCA of Aitchison distances ####
# pc.atdist=prcomp(atdist_matrix)

load('/home/baldanzi/Datasets/sleep_SCAPIS/pc_aitdist')

# Plotting the first 2 principal components 
# First: extracting the principal component value for every sample 
dat.plot=data.frame(pc.atdist$x, OSAcat = count$OSAcat)
dat.plot$OSAcat = factor(dat.plot$OSAcat, levels=c("no OSA", "Mild", "Moderate", "Severe"))

# Second: creating the scatter plot 
p1=ggplot(dat.plot,aes(x=PC1,y=PC2, color=OSAcat))+
  geom_point(size=1.1)+
  stat_ellipse(type = "t", size=1.3) +
  ggtitle("Aitchison distance") +
  xlab(paste0("PC1 \n (",round(summary(pc.atdist)$importance[2,1]*100,2),"% )")) +
  ylab(paste0("PC2 \n (",round(summary(pc.atdist)$importance[2,2]*100,2),"% )")) +
  theme(axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))

ggsave("Atdist.OSAcat.png", plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# Creating a Aitchison distance heatmap ####

paletteLength <- 500
myColor <- colorRampPalette(c("white","indianred3"))(paletteLength)

#annotation = data.frame(OSAcat = count$OSAcat)
annotation = data.frame(OSAcat = atdist_df$OSAcat)
rownames(annotation) = colnames(atdist_matrix)
annotation$OSAcat = factor(annotation$OSAcat, levels=c("no OSA", "Mild", "Moderate", "Severe"))

# Heatmap
p2=pheatmap(atdist_matrix,main="Heatmap",
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
atdist_df$OSAcat = factor(atdist_df$OSAcat, levels=c("no OSA","Mild","Moderate","Severe"))

# Box plot - Ait. Dist. between no OSA and the other categories
a = c(atdist_df[OSAcat=="no OSA",SCAPISid],"SCAPISid","OSAcat")
ADnoOSA= atdist_df[,a,with=F]
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
setDT(atdist_df)
a = c(atdist_df[OSAcat=="Severe",SCAPISid],"SCAPISid","OSAcat")
ADsevere= atdist_df[,a,with=F]
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

#-----------------------------------------------------------------------------#
# Relation to T90 (percentage of time with oxygen saturation < 90%) ####

valid.odi = fread("/home/baldanzi/Datasets/sleep_SCAPIS/validodi.MGS.Upp.tsv", sep = "\t")
valid.odi = valid.odi[!is.na(sat90),]

valid.odi[sat90<10,sat90cat :=1]
valid.odi[sat90>=10 & sat90<20,sat90cat :=2]
valid.odi[sat90>=20,sat90cat :=3]
valid.odi[,sat90cat:= factor(valid.odi$sat90cat, levels = c(1,2,3),
                             labels = c("<10","10-20",">20"))]
valid.odi %>% group_by(sat90cat) %>% summarise(N=length(sat90), mean=mean(sat90), median=median(sat90))
# A tibble: 3 x 4
#sat90cat       N         mean      median
#<fct>        <int>       <dbl>       <dbl>
#  1 <10       2558       1.90         1
#2 10-20       428        14.0        14
#3 >20         635        43.4        37

# Shannon diversity 
#T90 - Shannon diversity ####
a = grep("____",names(valid.odi),value=T) # vector with MGS names 
valid.odi$shannon=diversity(valid.odi[,a, with=F],index="shannon") #estimating shannon per individual

valid.odi %>% group_by(sat90cat) %>% summarise(N=length(shannon), mean=mean(shannon), median=median(shannon))

# Scatter plot: Shannon index against T90
a = with(valid.odi, cor(sat90, shannon, method = "spearman"))
a = round(a,3)
p1 = ggplot(data=valid.odi,aes(x=sat90, y= shannon)) + 
  geom_point() + ggtitle(paste0("Shannon index and T90 \nSpear. cor",a)) + 
  xlab("T90") +
  ylab("Shannon index") +
  geom_smooth(method='lm', alpha=.8) +
  theme(plot.title = element_text(hjust=.5, size=14),
        axis.title = element_text(size=14))

ggsave("scatter.shannon.T90.png", plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")
#Box plot
res=with(valid.odi, kruskal.test(shannon ~ sat90cat))
p2 = ggplot(data=valid.odi) + 
  geom_boxplot(aes(x=sat90cat, y= shannon, fill=sat90cat)) + 
  ggtitle(paste0("Shannon index by T90 groups\nKruskal-Wallis: p.value=",
                 formatC(res$p.value,format= "e", digits=2)))+
  xlab("T90") +
  ylab("Shannon index") +
  theme(title = element_text(hjust=.5, size=14),
        axis.title = element_text(size=12))

ggsave("box.shannon.T90.png", plot = p2, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")



# Bray curtis
# calculating BC for the samples with valid T90 measurements 
a=grep("____",names(valid.odi),value=T)
MGS=valid.odi[,a, with=F]
BC=as.matrix(vegdist(MGS,method="bray"))
rownames(BC)=colnames(BC)=valid.odi$SCAPISid

# PCA of BC
pc.bray=prcomp(BC) 

#  Extracting the proportion of variance explained for each principal component (eigenvalues)
prop=data.frame("PC"=names(summary(pc.bray)$importance[2,]),
                "value"=summary(pc.bray)$importance[2,]*100,
                "cumulative_prop" = summary(pc.bray)$importance[3,]*100,
                stringsAsFactors = T)

# Saving the BC matrix  ####
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
fwrite(BC, file = "T90.BCmatrix.csv", sep = ",")

# PCA of BC Biplot 
# First: extracting the principal component value for every sample 
dat.plot=data.frame(pc.bray$x, sat90cat = valid.odi$sat90cat)

# Second: creating the scatter plot ####
p4=ggplot(dat.plot,aes(x=PC1,y=PC2, color=sat90cat))+
  geom_point(size=1.1)+
  stat_ellipse(type = "t", size=1.3) +
  ggtitle("Bray-curtis dissimilarity - T90") +
  xlab(paste0("PC1 \n (",round(summary(pc.bray)$importance[2,1]*100,2),"% )")) +
  ylab(paste0("PC2 \n (",round(summary(pc.bray)$importance[2,2]*100,2),"% )")) +
  theme(axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))

ggsave("BC.T90.png", plot = p4, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")


# T90 - BC heatmap ####

paletteLength <- 500
myColor <- colorRampPalette(c("white","indianred3"))(paletteLength)
annotation = data.frame(sat90cat = valid.odi$sat90cat)
rownames(annotation) = colnames(BC)
annotation$sat90cat = factor(annotation$sat90cat, levels = c("<10","10-20",">20"))

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


# T 90 - Aitchison Distance ####

# Aitchison distance the Eucledian distance between CLR transformed count data. 

# Importing count data:
count = readRDS('/home/baldanzi/Datasets/MGS/upugut03.mgsCounts.rds')

# Zeros are replaced with 1 (there were no previous 1 in the data)
count[count==0]=1

#To fix some added characters to the names
count=as.data.frame(count)
rownames(count)[grep("gut.+_b",rownames(count),value=F)]=gsub("_b","",grep("gut.+_b",rownames(count),value=T))
rownames(count)[grep("gut.+b",rownames(count),value=F)]=gsub("b","",grep("gut.+b",rownames(count),value=T))
rownames(count)[grep("gut.+c",rownames(count),value=F)]=gsub("c","",grep("gut.+c",rownames(count),value=T))

# Adjusting the ID variable 
count$sample.id=rownames(count)

#Merging with the SCAPISid
selected_id = fread("/home/baldanzi/Datasets/scapis_idkey/selected_participants_id_sample_id.tsv", 
                    header=T, na.strings=c("", "NA"))
count = merge(count, selected_id, by="sample.id",all=F,all.x=T, all.y=F)

ids=fread("/home/baldanzi/Datasets/scapis_idkey/id_conversion.txt", header=T) 
ids$id = as.character(ids$export_id)
setnames(ids, "subject_id", "SCAPISid")
count = merge(count, ids, by = "id", all=F, all.x=T,all.y=F)

# Only keeping the individuals that have valid T90
count=count[count$SCAPISid %in% valid.odi[,SCAPISid],]

# Merging counts with valid.odi 
count = count %>% arrange(SCAPISid)
a = c("SCAPISid", "sat90cat", "sat90", "age", "Sex")
count=merge(count,valid.odi[,a,with=F], by="SCAPISid",all=F, all.x=F, all.y=T)

# Aitchison distance ####
rownames(count)=count$SCAPISid
mgs_count = as.matrix(count[,grep("HG3A",names(count))])

print("Aitchison distance calculation")
t0 = Sys.time()
atdist=aDist(mgs_count)
t1 = Sys.time()
print(paste0("time=" ,t1-t0))

# Transform aitchison distances into a matrix
atdist_matrix=as.matrix(atdist)

# Saving matrix 
fwrite(atdist_matrix, file = "T90.ADmatrix.csv", sep = ",")

# PCA of Aitchison distances ####
pc.atdist=prcomp(atdist_matrix)
# Plotting the first 2 principal components 
# First: extracting the principal component value for every sample 
dat.plot=data.frame(pc.atdist$x, sat90cat = count$sat90cat)

# Second: creating the scatter plot 
p1=ggplot(dat.plot,aes(x=PC1,y=PC2, color=sat90cat))+
  geom_point(size=1.1)+
  stat_ellipse(type = "t", size=1.3) +
  ggtitle("T90 - Aitchison distance") +
  xlab(paste0("PC1 \n (",round(summary(pc.atdist)$importance[2,1]*100,2),"% )")) +
  ylab(paste0("PC2 \n (",round(summary(pc.atdist)$importance[2,2]*100,2),"% )")) +
  theme(axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))

ggsave("Atdist.T90.png", plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# T90 - Ait Dist Heatmap
# T90 - BC heatmap ####

paletteLength <- 500
myColor <- colorRampPalette(c("white","indianred3"))(paletteLength)
annotation = data.frame(sat90cat = valid.odi$sat90cat)
rownames(annotation) = colnames(atdist_matrix)
annotation$sat90cat = factor(annotation$sat90cat, levels = c("<10","10-20",">20"))
rownames(atdist_matrix) = colnames(atdist_matrix)

# Heatmap
p1= pheatmap(atdist_matrix,main="Heatmap",
             legend_labels = "Aitchison Distance - T90",
             color = myColor,
             show_rownames = FALSE, show_colnames = FALSE,
             cluster_rows = T,cluster_cols = T,
             annotation_row = annotation,
             annotation_col = annotation)

ggsave("T90.AD.heatmap.png", plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")
#----------------------------------------------------------------------------#

# Comparison of relative abundance vs CLR transformed counts ####

#### REMOVE THIS ####
valid.ahi= fread('validsleep_MGS.shannon.BC_Upp.tsv', head=T)
load("data_table1")
a = names(dat1[,-"SCAPISid"])
valid.ahi[,(a):=NULL]
valid.ahi = merge(dat1,valid.ahi, by = "SCAPISid", all=T, all.x=F, all.y=F)
count.clr = fread("/home/baldanzi/Datasets/sleep_SCAPIS/validsleep.mgs_clr_count.tsv",sep='\t')
#---###---###---###---###

# Extracting mgs relative abundances 
a = grep("___",names(valid.ahi),value=T)
mgs.ra=as.matrix(valid.ahi[,a,with=F])
rownames(mgs.ra) = valid.ahi[,SCAPISid]

# Extracting mgs clr 
a = grep("HG3A",names(count.clr),value=T)
mgs.clr = as.matrix(count.clr[,a,with=F])
rownames(mgs.clr) = count.clr[,SCAPISid]

# Correlation between the two mgs transformation by individuals 
cor.ind=sapply(1:nrow(mgs.clr), function(i){
  cor(mgs.clr[i,],mgs.ra[i,], method = "spearman")})

summary(cor.ind)
# Min.    1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.9998  0.9999  1.0000   1.0000  1.0000   1.0000 

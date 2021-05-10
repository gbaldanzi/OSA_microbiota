# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics continued : script to DESCRIBE alpha and beta-diversity

# This script will calculate alpha diversith and create plots and table describing alpha
#diversity in relation AHI and OSA, as well as other covarites. 

# Descriptive Statistics 

# Loading packages
library(tidyverse)
library(Hmisc)
library(compareGroups)
library(flextable)
library(data.table)
library(vegan)
library(pheatmap)
library(RColorBrewer)
library(robCompositions)

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

summary(valid.ahi[,shannon])
valid.ahi[,mean(shannon), by=OSAcat]
valid.ahi[,median(shannon), by=OSAcat]

# Compare the Shannon index across two groups using Kruskal-Wallis
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
BC=as.matrix(vegdist(MGS,method="bray"))
rownames(BC)=colnames(BC)=valid.ahi$SCAPISid

# Principal coordinates analysis based on BC 
pc.bray=prcomp(BC)

# Extracting the proportion of variance explained for each principal component (eigenvalues)
prop=data.frame("PC"=names(summary(pc.bray)$importance[2,]),
                "value"=summary(pc.bray)$importance[2,]*100,
                "cumulative_prop" = summary(pc.bray)$importance[3,]*100,
                stringsAsFactors = T)

# Saving data with shannon index and principal components of BC 
# Saving datasets ####
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
save(pc.bray, file = 'pc_BC')
fwrite(valid.ahi, file = "validsleep_MGS.shannon.BC_Upp.tsv", sep = "\t")
fwrite(BC, file = "OSA.BCmatrix.csv", sep = ",")

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

# Creating a BC heatmap ####

paletteLength <- 500
myColor <- colorRampPalette(c("white","indianred3"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(0, 0.5, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))

annotation = data.frame(OSAcat = valid.ahi[,OSAcat])
rownames(annotation) = valid.ahi$SCAPISid

# Heatmap
p5=pheatmap(BC,main="Heatmap",
            legend_labels = "Bray-Curtis Dissimilarity",
            color = myColor,show_rownames = FALSE, show_colnames = FALSE,
            # breaks = myBreaks,
            cutree_rows = 2,
            cutree_cols = 2,cluster_rows = F,cluster_cols = F, 
            annotation_row = annotation,
            annotation_col = annotation)

ggsave("BC.heatmap", plot = p5, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

#-------------------------------------------------------------------------#
# Beta diversity based on Aitchison distance. ####

# Aitchison distance the Eucledian distance between CLR transformed count data. 

# Importing count data:
count = readRDS('/home/baldanzi/Datasets/MGS/upugut03.mgsCounts.rds')

# Summary of total counts per sample 
row_sum=rowSums(count)
summary(row_sum)
#  Min. 1st Qu.  Median    Mean  3rd Qu.    Max. 
# 1397  587162  657876   672899  740013   2016207 

# Number of zeros in the data
sum(count==0)
# 8142757

# Zeros are replaced with 1 (there were no preivous 1 in the data)
count[count==0]=1

# Central log-transformation of count data 
count.clr = cenLR(count)
count=as.data.frame(count)
count.clr = as.data.frame(count.clr$x.clr)


#To fix some added characters to the names
rownames(count)[grep("gut.+_b",rownames(count),value=F)]=gsub("_b","",grep("gut.+_b",rownames(count),value=T))
rownames(count)[grep("gut.+b",rownames(count),value=F)]=gsub("b","",grep("gut.+b",rownames(count),value=T))
rownames(count)[grep("gut.+c",rownames(count),value=F)]=gsub("c","",grep("gut.+c",rownames(count),value=T))

rownames(count.clr)[grep("gut.+_b",rownames(count.clr),value=F)]=gsub("_b","",grep("gut.+_b",rownames(count.clr),value=T))
rownames(count.clr)[grep("gut.+b",rownames(count.clr),value=F)]=gsub("b","",grep("gut.+b",rownames(count.clr),value=T))
rownames(count.clr)[grep("gut.+c",rownames(count.clr),value=F)]=gsub("c","",grep("gut.+c",rownames(count.clr),value=T))


# Adjusting the ID variable 
count$sample.id=rownames(count)
count.clr$sample.id=rownames(count.clr)



#Merging with the SCAPISid
selected_id = fread("/home/baldanzi/Datasets/scapis_idkey/selected_participants_id_sample_id.tsv", 
                    header=T, na.strings=c("", "NA"))
count = merge(count, selected_id, by="sample.id",all=F,all.x=T, all.y=F)
count.clr = merge(count.clr, selected_id, by="sample.id",all=T,all.x=F, all.y=F)

ids=fread("/home/baldanzi/Datasets/scapis_idkey/id_conversion.txt", header=T) 
ids$id = as.character(ids$export_id)
setnames(ids, "subject_id", "SCAPISid")
count = merge(count, ids, by = "id", all=F, all.x=T,all.y=F)
count.clr = merge(count.clr, ids, by = "id", all=F, all.x=T,all.y=F)

# Only keeping the individuals that have valid sleep recording
count=count[count$SCAPISid %in% valid.ahi[,SCAPISid],]

# Merging counts with valid.ahi 
a = c("SCAPISid", "OSAcat", "ahi", "age", "Sex")
count=merge(count,valid.ahi[,a,with=F], by="SCAPISid",all=F, all.x=F, all.y=T)

# Aitchison distance ####
rownames(count)=count$SCAPISid
mgs_count = as.matrix(count[,grep("HG3A",names(count))])


t0 = Sys.time()
atdist=aDist(mgs_count)
t1 = Sys.time()
print(paste0("time=" ,t1-t0))

# Transform aitchison distances into a matrix
atdist_matrix=as.matrix(atdist)

# PCA of Aitchison distances ####
pc.atdist=prcomp(atdist_matrix)

# Plotting the first 2 principal components 
# First: extracting the principal component value for every sample 
dat.plot=data.frame(pc.atdist$x, OSAcat = count$OSAcat)

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
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(0, 0.5, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))

annotation = data.frame(OSAcat = count$OSAcat)
rownames(annotation) = count$SCAPISid

# Heatmap
p2=pheatmap(atdist_matrix,main="Heatmap",
            legend_labels = "Aitchison Distance",
            color = myColor,show_rownames = FALSE, show_colnames = FALSE,
            # breaks = myBreaks,
            cutree_rows = 2,
            cutree_cols = 2,cluster_rows = F,cluster_cols = F, 
            annotation_row = annotation,
            annotation_col = annotation)

ggsave("AitDist.heatmap.png", plot = p2, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# Merge clr-data with valid sleep recordings 
a = c("SCAPISid", "OSAcat", "ahi", "age", "Sex")
count.clr=merge(count.clr,valid.ahi[,a,with=F], by="SCAPISid",all=F, all.x=F, all.y=T)

# Saving clr-transformed data 
fwrite(count.clr,"/home/baldanzi/Datasets/sleep_SCAPIS/validsleep.mgs_clr_count.tsv",sep='\t')

#Saving the atchison distances 
atdist_df = as.data.frame(atdist)
atdist_df$SCAPISid = count$SCAPISid
atdist_df$OSAcat = count$OSAcat
fwrite(atdist_df,"/home/baldanzi/Datasets/sleep_SCAPIS/OSA.aitchison_distmatrix.csv",sep=',')

# Saving the principal components analysis 
save(pc.atdist, file = '/home/baldanzi/Datasets/sleep_SCAPIS/pc_aitdist')

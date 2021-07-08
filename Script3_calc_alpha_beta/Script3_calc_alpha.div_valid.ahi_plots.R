# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will create plots and table describing alpha in participants with valid ahi

# Loading packages
library(grid)

  setwd("/home/baldanzi/Datasets/sleep_SCAPIS")

# Uploading dataset with variables of class "factor"
  load("data_table1")
  
  
# Importing data on individuals with valid.ahi measurement 
  valid.ahi <- fread("validsleep_MGS.shannon_Upp.tsv", sep = "\t")
  
  
#merging with the data set that has variables as factors 
  a = names(dat1[,-"SCAPISid"])
  valid.ahi[,(a):=NULL]
  valid.ahi = merge(dat1,valid.ahi, by = "SCAPISid", all=T, all.x=F, all.y=F)


# Compare the Shannon index across OSA groups using Kruskal-Wallis #### 
res=with(valid.ahi, kruskal.test(shannon ~ OSAcat))

# Histogram of Shannon diversity 
  p3_1=ggplot(data=valid.ahi,aes(x=shannon))+
        geom_histogram()+
        ggtitle("Histogram - Shannon diversity")
  
  ggsave("hist.shannon", plot = p3_1, device = "png", 
         path = "/home/baldanzi/Sleep_apnea/Descriptive/")
  
# Scatter plot: Shannon index against AHI 
  a = cor(valid.ahi$ahi, valid.ahi$shannon, method = "spearman")
  
  p1 = ggplot(data=valid.ahi,aes(x=ahi, y= shannon)) + 
    geom_point(color="lightskyblue2") + 
    ggtitle(paste0("Shannon index and AHI \n Spearman cor.= ",round(a,3))) + 
    xlab("Apnea-hypopnea index") +
    ylab("Shannon index") +
    geom_smooth(method='lm', alpha=.8) +
    geom_vline(xintercept = 5, color = "red", linetype = "twodash") + 
    geom_vline(xintercept = 15, color = "red", linetype = "twodash") + 
    geom_vline(xintercept = 30, color = "red", linetype = "twodash") +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  ggsave("scatter.shannon.ahi.png", plot = p1, device = "png", 
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

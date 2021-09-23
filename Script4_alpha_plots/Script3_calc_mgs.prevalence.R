# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will calculate MGS prevalences  

# Descriptive Statistics 

rm(list=ls())

# Loading packages
library(Hmisc)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")

# Import pheno data
valid.ahi <- readRDS("validsleep_MGS.shannon_Upp.rds")

# MGS prevalence ####

# MGS prevalence is the number of participants in which a certain MGS is present

# presence-absence transformation:
# If a species is present, it is transformed to 1
# If it is absent, it is transformed to zero
mgs = grep("____H",names(valid.ahi))
data_pa <- decostand(x = valid.ahi[,mgs,with=F], "pa")

# calculate sum per species
  data_sum <- apply(data_pa, 2, sum)

# Histogram of MGS prevalence
  p1 = ggplot(data = as.data.frame(data_sum), aes(x=data_sum)) + 
              geom_histogram(bins = 70,color="black", fill="lightskyblue2") +
              geom_density(alpha=.3, fill="#FF6666") +
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
              theme_light()+
              theme(plot.title = element_text(hjust = .5,size = 14,face = "bold"),
                    axis.title.x = element_text(size=14))

  ggsave("hist.MGS.prevalence",plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")


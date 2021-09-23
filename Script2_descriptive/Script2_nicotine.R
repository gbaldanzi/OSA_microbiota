# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# Last update: 2021-09-22

library(data.table)
library(tidyverse)

  # Import data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  valid.ahi <- pheno[valid.ahi=='yes',]

# Because there was a high number of individuals with missing information on 
#smoking status, we decided to look more into Nicotine metabolites from metabolomics
#measurements as a proxy for smoking status. 

# Nicotine metabolites ####
# Importing data on nicotine metabolites
  nicotine = fread("/home/baldanzi/Datasets/Metabolon/nicotine_metab.csv", header = T, 
                 na.strings = c("NA",""))
  nicot_label = c("cotinine", "hydroxycotinine", "cotinine N-oxide")
  label(nicotine$MET_848) = nicot_label[1]
  label(nicotine$MET_100002717) = nicot_label[2]
  label(nicotine$MET_100002719) = nicot_label[3]

# How many individuals have NA values for each of the nicotine metabolites 
  apply(nicotine, 2, function(x) {sum(is.na(x))}) # 28,28,95


#box plots of nicotine metabolites by smokestatus 
#Cotinine 
nic1 = ggplot(nicotine, aes(y = MET_848,x = smokestatus, fill= smokestatus)) + 
  geom_boxplot() + ggtitle("Cotinine and smoking status") + ylab("Cotinine") +
  theme(title = element_text(hjust = 0.5, size=16, face = "bold"), 
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=13))
ggsave("smokcotinine.png", plot = nic1, device = "png", path="/home/baldanzi/Sleep_apnea/Descriptive/")

#Hydroxycotinine 
nic2 = ggplot(nicotine, aes(y = MET_100002717,x = smokestatus, fill= smokestatus)) + 
  geom_boxplot() + ggtitle("Hydroxycotinine and smoking status") + ylab("Hydroxycotinine") +
  theme(title = element_text(hjust = 0.5, size=16, face = "bold"), 
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=13))
ggsave("smokhydroxycotinine.png", plot = nic2, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")

# Cotinine N-oxide
nic3 = ggplot(nicotine, aes(y = MET_100002719,x = smokestatus, fill= smokestatus)) + 
  geom_boxplot() + ggtitle("Cotinine N-oxide and smoking status") + ylab("Cotinine N-oxide") +
  theme(title = element_text(hjust = 0.5, size=16, face = "bold"), 
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=13))
ggsave("smokcotinineNoxide.png", plot = nic3, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")

# Histograms for Cotinine by smoking status 
nic4 = ggplot(data = nicotine, aes(x=MET_848)) +
  geom_histogram(data=subset(nicotine,smokestatus=="never"), fill="purple", alpha=0.6)+
  ggtitle("Never smokers") + xlab("Cotinine")
nic5 = ggplot(data = nicotine, aes(x=MET_848)) +
  geom_histogram(data=subset(nicotine,smokestatus=="former"), fill="blue", alpha=0.6)+
  ggtitle("Former smokers") + xlab("Cotinine")
nic6 = ggplot(data = nicotine, aes(x=MET_848)) +
  geom_histogram(data=subset(nicotine,smokestatus=="current"), fill="green", alpha=0.6)+
  ggtitle("Current smokers") + xlab("Cotinine")
nic456 = ggarrange(nic4, nic5, nic6, nrow=1)

ggsave("smokcotinine_hist.png", plot = nic456, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")

# Histogram for Cotinine by smoking status zoomed at the origin  
nic4 = ggplot(data = nicotine[MET_848<1,], aes(x=MET_848)) +
  geom_histogram(data=subset(nicotine[MET_848<1,],smokestatus=="never"), fill="purple", alpha=0.6)+
  ggtitle("Never smokers") + xlab("Cotinine")
nic5 = ggplot(data = nicotine[MET_848<1,], aes(x=MET_848)) +
  geom_histogram(data=subset(nicotine[MET_848<1,],smokestatus=="former"), fill="blue", alpha=0.6)+
  ggtitle("Former smokers") + xlab("Cotinine")
nic6 = ggplot(data = nicotine[MET_848<1,], aes(x=MET_848)) +
  geom_histogram(data=subset(nicotine[MET_848<1,],smokestatus=="current"), fill="green", alpha=0.6)+
  ggtitle("Current smokers") + xlab("Cotinine")
nic456 = ggarrange(nic4, nic5, nic6, nrow=1)

ggsave("smokcotinine_hist_zoomed.png", plot = nic456, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")

# Cross tabulation Cotinine using min value
# Cotinine was transformed into a factor variables (presente or absent)
#depending if it was above or below the min value detected. 
nicotine[MET_848>0.0003,cotinine_pa:='present']
nicotine[MET_848<=0.0003,cotinine_pa:='absent']

t = ctable(nicotine$cotinine_pa, nicotine$smokestatus, round.digits = 0, prop = "c")

# ROC cotinine and current smoking 
# Trying to find the best Cotinine values to distinguish smokers from non-smokers. 
library(pROC)
nicotine[smokestatus=="current", current :=1]
nicotine[smokestatus!="current", current :=0]
my_roc <- roc(nicotine$current, nicotine$MET_848)
a = coords(my_roc, "best", best.method = "youden", 
           ret = c("threshold", "specificity", "sensitivity"), 
           transpose=TRUE)
a
png('/home/baldanzi/Sleep_apnea/Descriptive/cotinineROC.png', width=600, height = 600)
par(cex.axis=1.3, cex.main=1.4 )
plot(my_roc) 
title(paste0("ROC-smoking and cotinine\nAUC=",round(auc(my_roc),3)))
dev.off()

# Cross tabulation Cotinine using ROC value
# Cotinine was transformed into a factor variables (presente or absent)
#depending if it was above or below the threshold value determined by the ROC
nicotine[,cotinine_pa:=NULL]
nicotine[MET_848>a[1],cotinine_pa:='present']
nicotine[MET_848<=a[1],cotinine_pa:='absent']
t = with(nicotine, table(cotinine_pa, smokestatus))
t
round(prop.table(t, 2),3)*100

# Mean and median cotinine levels by smoking status and number of individuals with cotinine
#level above ROC threshold. 
nicotine[!is.na(MET_848),] %>%  group_by(smokestatus) %>% 
  summarise(N= length(MET_848), mean=mean(MET_848, na.rm = T),median = median(as.numeric(MET_848), na.rm=T), 
            cotinine_above_cutoff = sum(MET_848>a[1]))


# Snus use and cotinine serum level ####
n = nrow(nicotine)
snus1 = ggplot(nicotine, aes(y = MET_848,x = snus_ever1mo, fill= snus_ever1mo)) + 
  geom_boxplot() + ggtitle(paste0("Cotinine and snus\nn=",n)) + ylab("Cotinine") +
  theme(title = element_text(hjust = 0.5, size=14, face = "bold"), 
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=13))
ggsave("snuscotinine.png", plot = snus1, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")

nicotine[!is.na(MET_848),] %>%  group_by(snus_ever1mo) %>% 
  summarise(N= length(MET_848), mean=mean(MET_848, na.rm = T),
            median = median(as.numeric(MET_848), na.rm=T))

# Relation between snus use and cotinine serum level excluding smokers
n = nrow(nicotine[smokestatus!="current",])
snus2 = ggplot(nicotine[smokestatus!="current",], 
               aes(y = MET_848,x = snus_ever1mo, fill= snus_ever1mo)) + 
  geom_boxplot() + ggtitle(paste0("Cotinine and snus in never smokers\nn=",n)) + ylab("Cotinine") +
  theme(title = element_text(hjust = 0.5, size=14, face = "bold"), 
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=13))
ggsave("snuscotinine_neversmokers.png", plot = snus2, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")

nicotine[!is.na(MET_848) & smokestatus!="current",] %>%  group_by(snus_ever1mo) %>% 
  summarise(N= length(MET_848), mean=mean(MET_848, na.rm = T),median = median(as.numeric(MET_848), na.rm=T), 
            above.threshold = sum(MET_848>a[1]))


# Smoking status based on MET_848 (Cotinine) excluding snus users. 
nicotine[!is.na(MET_848) & snus_ever1mo!="YES",] %>%  group_by(smokestatus) %>% 
  summarise(N= length(MET_848), mean=mean(MET_848, na.rm = T),median = median(as.numeric(MET_848), na.rm=T), 
            above.threshold = sum(MET_848>a[1]))


nrow(nicotine[SCAPISid %in% valid.ahi[!is.na(cqto015) & is.na(smokestatus),SCAPISid],])

# Missingness in snus use 
sum(is.na(valid.ahi[,cqto015])) # 128 

# Cross tabulation snus and smoking ####
snus_smoking = data.frame(snus_missing = numeric(length = nrow(valid.ahi)), 
                          smoke_missing = numeric(length = nrow(valid.ahi)), 
                          snus=valid.ahi$cqto015,
                          smoke=valid.ahi$smokestatus)
snus_smoking$snus_missing = ifelse(is.na(valid.ahi$cqto015), 0, 1)
snus_smoking$smoke_missing = ifelse(is.na(valid.ahi$smokestatus), 0, 1)
snus_smoking[,c("snus_missing","smoke_missing")] = apply(snus_smoking[,c("snus_missing","smoke_missing")]
                                                         , 2, function(x){
                                                           factor(x, levels=c(0,1), labels = c("NA", "not NA"))})
label(snus_smoking$snus_missing) = "Snus"
label(snus_smoking$smoke_missing) = "Smoking status"

save(snus_smoking, file='/home/baldanzi/Sleep_apnea/Descriptive/snus_smoking')
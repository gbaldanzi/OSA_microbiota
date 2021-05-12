# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

print("#Descriptive Statistics#")
print(Sys.time())

# Loading packages
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(compareGroups)
library(flextable)
library(data.table)
library(vegan)
library(pheatmap)
library(RColorBrewer)
library(rio)

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")

# Uploading dataset with phenotype variables for initial descriptive statistics 
load("data_table1")

#Variables are: OSAcat , age, Sex, smokestatus, Alkohol, BMI, WaistHip, educat,
# leisurePA, placebirth, diabd, hypertension, dyslipidemia, diabmed, 
# hypermed, dyslipmed, PPI, Fibrer, ahi, odi, sat90, ESS, cpap, splint

# Variable diagnostics ####
# Number of individuals 
nrow(dat1) #3206

# Checking for missingness in the variables. 
miss = data.frame(Missing = apply(dat1,2,function(x){sum(is.na(x))}))
a = apply(dat1,2,function(x){sum(is.na(x))/length(x)})
miss = cbind(Variables = rownames(miss), miss, Percentage = paste0(round(a*100,0)," %"))
flextable(miss)

# Number of individuals having at least one missing value for the covariates
perrow = apply(dat1, 1, function(x) {any(is.na(x))})
sum(perrow) 

# Summary of variables
summary = sapply(names(dat1), function(x) {summary(dat1[[x]])})

# "Table 1" - Population characteristics ####
datatable1 = as.data.frame(dat1[,-1])  # Passing data to a new object without the SCAPIS ID

# For the purpose of this table, individuals with missing values for the following variables
# were included in the "no" category. (only category "yes" is displayed in the table)
for(i in which(names(datatable1) %in% c('hypertension', 'dyslipidemia', 'diabmed',
                                'hypermed','dyslipmed', 'ppi','cpap', 'splint'))){
 datatable1[is.na(datatable1[i]),i] = "no"
}

#Creating labels for the variables 
mylabel <- c("OSA severity", "Age (yrs)", "Sex", "Smoking status", "Alcohol intake (g)", 
             "BMI (kg/m2)", "WHR" ,"Highest education achieved", "Self-reported leisure physical activity", 
             "Place of birth", "Type 2 diabetes", "Hypertension", "Dyslipidemia",
             "Diabetes medication",
             "Hypertension medication", "Dyslipidemia medication", "PPI", "Fiber", "AHI (events/h)", "ODI (events/h)", 
             "T90%", "ESS", "CPAP", "Oral appliance")

j=1
for(i in names(datatable1)){
  label(datatable1[[i]]) <- mylabel[j]
  j=j+1
}

# Create table of population characteristics by OSA group (no OSA, Mild, Moderate, or Severe OSA)
t = compareGroups(OSAcat ~ age + Sex + smokestatus + Alkohol + BMI +WaistHip+ educat +
                  leisurePA + pob + diabd + hypertension + dyslipidemia + diabmed + 
                  hypermed + dyslipmed+ppi+Fibrer+ ahi+ odi+ sat90+ ESS + cpap+ splint,
                  data= datatable1, 
                  include.miss = TRUE, chisq.test.perm = TRUE)
t1 = createTable(t, hide.no = "no")
t1
export2html(t1, file='/home/baldanzi/Sleep_apnea/Descriptive/sleepapnea_table1.html', header.labels = c(p.overall = "p-value"))

rm(datatable1)

#______________________________________________________________________________

# Descriptive statistics ####

setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
# Importing data on individuals with valid.ahi measurement (phenotype+MGS)
valid.ahi = fread("validsleep.MGS.Upp.tsv", header=T, na.strings=c("", "NA"))

#merging with the data set that have variables as factor variables
a = names(dat1[,-"SCAPISid"]) # Assigning to object the column names from dat1
valid.ahi[,(a):=NULL] # Erasing from the main dataset the columns to be substituted with factor variables 
valid.ahi = merge(dat1,valid.ahi, by = "SCAPISid", all=T, all.x=F, all.y=F) 


#-----------------------------------------------------------------------#
# Because there was a high number of individiausl with missing information on 
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

#---------------------------------------------------------------------------#
# We wanted to check if blood samples for metabolomics and Sleep records (ApneaLink)
#took place close to each other (difference in calender data between the two)


#### Dates ####
# Was sleep apnea assessed close to when metabolomics was collected? 
valid.ahi[,metabolon_collection_date:=as.Date(valid.ahi[,metabolon_collection_date], 
                                              format = "%d-%b-%y")]
a = c("date", "metabolon_collection_date", "AnthropometryCollectionDate","LdlSamplingDate" )
# Number of participants for which metabolon_collection_date is missing 
nrow(valid.ahi[is.na(metabolon_collection_date),]) # 1246

# Number of participants for which metabolon_colletion and sleep apnea assessment date 
# are the same
nrow(valid.ahi[date==metabolon_collection_date,]) # 1511

# range and mean difference between the two dates 
diff=abs(as.Date(valid.ahi[,date])-valid.ahi[,metabolon_collection_date])
summary(as.numeric(diff))
nrow(valid.ahi[diff>30,a,with=F]) # 10 individuals have a difference greater than 
# 30 days between dates

# Sleep assessment dates and Antropometric Collection Date 
# Missing data in anthropometric collection 
sum(is.na(valid.ahi[,AnthropometryCollectionDate])) # ZERO
# difference between the two dates 
diff2=abs(as.Date(valid.ahi[,date])-as.Date(valid.ahi[,AnthropometryCollectionDate]))
summary(as.numeric(diff2))
nrow(valid.ahi[diff2>30,a,with=F]) # 23 had a difference greater than 30 days

#--------------------------------------------------------------------------#
#### Missingness #### 
# Variables of interest 
listvar = c("educat", "leisurePA", "hypertension","smokestatus", 
            "diabd", "pob", "Alkohol", "ppi")

# Selecting 
contain_missing = which(apply(dat1[,listvar,with=F], 1, function(x){any(is.na(x))}))
dat1[,incomplete_obs:=NULL]
dat1[contain_missing,incomplete_obs:=1]
dat1[!contain_missing,incomplete_obs:=2]
dat1$incomplete_obs = factor(dat1[,incomplete_obs], 
                               levels= c(1,2),
                               labels = c("Incomplete obs.", 
                                          "Complete obs."))

# Preparing data for table: Participants characteritics divided by 
#observation complete or incomplete 
datafortable = as.data.frame(dat1[,-1]) 

# Labelling variables 
mylabel <- c("OSA severity", "Age (yrs)", "Sex", "Smoking status", "Alcohol intake (g)", 
             "BMI (kg/m2)", "Highest education achieved", "Self-reported leisure physical activity", 
             "Place of birth", "Type 2 diabetes", "Hypertension", "Dyslipidemia",
             "Diabetes medication",
             "Hypertension medication", "Dyslipidemia medication", "Fiber", "AHI (events/h)", "ODI (events/h)", 
             "%time Sat=<90%", "CPAP", "Oral appliance", "Incomplete observation")
j=1
for(i in names(datafortable)){
  label(datafortable[[i]]) <- mylabel[j]
  j=j+1
}

# Creating the table 
t = compareGroups(incomplete_obs ~ ahi + age + Sex + BMI+  Alkohol + smokestatus + educat +
                    leisurePA + pob + diabd + hypertension, data= datafortable, 
                  include.miss = TRUE, chisq.test.perm = TRUE)
t1 = createTable(t)
t1
export2html(t1, file='/home/baldanzi/Sleep_apnea/Descriptive/completVSincompletObs.html', 
            header.labels = c(p.overall = "p-value"))
save(t1, file='/home/baldanzi/Sleep_apnea/Descriptive/completVSincompletObs')

#---------------------------------------------------------------------------#
#### Clinical microbiomics variables ####
# Variables related to how fecal samples were analysed to produce the gut microbiota data

#Number of plates 
length(unique(valid.ahi$plate)) # 55 plates

# Number of occasions that the data was received
length(unique(valid.ahi$received))

# Average read.pairs.per.sample + dna.yield
n3 = mean(valid.ahi[,read.pairs.per.sample], na.rm = T)
n4 = mean(valid.ahi[,dna.yield], na.rm = T)

# 
tcm = compareGroups(OSAcat~read.pairs.per.sample + dna.yield, data= valid.ahi)
t1 = createTable(tcm)

#---------------------------------------------------------------------------#
#### Taxonomy ####
# Describe the taxonomic composition by different OSA severity groups. 

# Importing data with taxonomic information for every MGS 
taxonomy = import('/home/baldanzi/Datasets/MGS/MGS_taxonomic_information.tsv')

# Transforming taxonomic levels into factor variables
a = c("species", "genus", "family", 'order', 'class', 'phylum')
for(i in a){
taxonomy[[i]]=as.factor(taxonomy[[i]])
}
taxonomy$mgs2 = paste0(taxonomy$maintax,"____",taxonomy$mgs)

# Subsetting the data to create a data.frame containing only SCAPISid, ahi, OSAcat
#(OSA severity category), and the relative abundance by phylum for each participants 
a= c("SCAPISid","ahi","OSAcat",levels(taxonomy$phylum))
phylum_abundance=matrix(nrow = nrow(valid.ahi), ncol = length(a))
colnames(phylum_abundance)=a
phylum_abundance=as.data.frame(phylum_abundance)

a=c("SCAPISid","ahi","OSAcat")
phylum_abundance[,a]=as.data.frame(valid.ahi[,a,with=F])

for(i in levels(taxonomy$phylum)){ 
  a = taxonomy$mgs2[taxonomy$phylum==i]
    phylum_abundance[[i]]=rowSums(valid.ahi[,a,with=F])
}

# Determining the mean relative abundance of phylum by OSAcat
phylum_mean= phylum_abundance %>% gather(phyla_name,abundance, levels(taxonomy$phylum)) 
phylum_mean= phylum_mean %>% group_by(OSAcat,phyla_name) %>% summarise(average=mean(abundance))
phylum_mean= as.data.frame(phylum_mean)

#Transforming the OSAcat variable into a factor variable (important for the plot)
phylum_mean$OSAcat = factor(phylum_mean$OSAcat, levels=c("no OSA", "Mild", "Moderate", "Severe"))

# There are 14 phyla. To many for a stacked bar plot. 
# Phyla that have an average below the p25 for all categories will be transformed to "others"
a = quantile(phylum_mean$average,c(.25))
for(i in unique(phylum_mean$phyla_name)){
if(all(phylum_mean$average[phylum_mean$phyla_name==i]<a)){
  phylum_mean$phyla_name[phylum_mean$phyla_name==i]="other"
}
}
# Creating the bar plot stacked for phyla relative abundance by OSA cat 
mycolors <- colorRampPalette(brewer.pal(12, "Set2"))(14)

p1 = phylum_mean %>% ggplot(aes(x=OSAcat, y=average, fill=phyla_name)) +
  geom_bar(position = "stack", stat = "identity", color="black",lwd=.2) +
  ggtitle("Phyla abundance by sleep apnea severity") +
  ylab("Relative abudance") + xlab("Sleep apnea severity") + 
  theme(title=element_text(hjust=.5, size=12),
        axis.title = element_text(size=12)) +
  scale_fill_manual(name= "phyla", values = mycolors)

#Saving the plot
ggsave("phyla_abundance_OSAcat.png", plot = p1, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")

# Subsetting the data to create a data.frame containing only SCAPISid, ahi, OSAcat
#(OSA severity category), and the relative abundance by family for each participants 
a= c("SCAPISid","ahi","OSAcat",levels(taxonomy$family))
family_abundance=matrix(nrow = nrow(valid.ahi), ncol = length(a))
colnames(family_abundance)=a
family_abundance=as.data.frame(family_abundance)

a=c("SCAPISid","ahi","OSAcat")
family_abundance[,a]=as.data.frame(valid.ahi[,a,with=F])

for(i in levels(taxonomy$family)){ 
  a = taxonomy$mgs2[taxonomy$family==i]
    family_abundance[[i]]=rowSums(valid.ahi[,a,with=F])
}

# Determining the mean relative abundance of family by OSAcat
family_mean = family_abundance %>% gather(family_name,abundance, levels(taxonomy$family)) %>%
  group_by(OSAcat,family_name) %>% summarise(average=mean(abundance))
family_mean = as.data.frame(family_mean)


#Transforming the OSAcat variable into a factor variable (important for the plot)
family_mean$OSAcat = factor(family_mean$OSAcat, levels=c("no OSA", "Mild", "Moderate", "Severe"))

# There are 58 phyla. To many for a stacked bar plot. 
# Famlies that have an average below the p60 for all categories will be transformed to "others"
a = quantile(family_mean$average,c(.70))
for(i in unique(family_mean$family_name)){
  if(all(family_mean$average[family_mean$family_name==i]<a)){
    family_mean$family_name[family_mean$family_name==i]="other"
  }
}

# Creating the bar plot stacked for phyla relative abundance by OSA cat 
mycolors <- colorRampPalette(brewer.pal(12, "Set2"))(20)

p2 = family_mean %>% ggplot(aes(x=OSAcat, y=average, fill=family_name)) +
  geom_bar(position = "stack", stat = "identity", color="black",lwd=.2) +
  ggtitle("Family abundance by sleep apnea severity") +
  ylab("Relative abudance") + xlab("Sleep apnea severity") + 
  scale_fill_discrete(name = "families") +
  theme(title=element_text(hjust=.5, size=12),
        axis.title = element_text(size=12)) +
  scale_fill_manual(name= "families", values = mycolors)

#Saving the plot
ggsave("family_abundance_OSAcat.png", plot = p2, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")
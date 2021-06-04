# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script for data management and merging the sleep and MGS datasets

# Data management ####

# Loading packages
pacman::p_load(tidyverse, grid, chron, rio, Hmisc, sjmisc, summarytools, data.table)

rm(list=ls())

# Import sleep recording data 

setwd("/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording")
sleep = fread("Data_Sleep_SCAPIS_Uppsala.csv", header=T, na.strings=c("", "NA"))


# Data summary 
attach(sleep)
nr.idnr = length(idnr) # total number individuals 
nr.idnr.valid = length(idnr[flutv4h==1 & o2utv4h==1]) # number of ind with at least 4h of flow and sat monitoring 
fl = summary(fldeutv_min) 
fl = substr(times((fl%/%60 +  fl%%60 /60)/24), 1, 5) # Mean and median duration of flow monitoring 

sat = summary(spo2utv_min)
sat = substr(times((sat%/%60 +  sat%%60 /60)/24), 1, 5) # Mean and median duration of flow monitoring
detach(sleep)

# Date of recording ####

# Date interval 
range(as.Date(sleep$date, format = "%Y-%m-%d"),na.rm=T) #highest value with year = 2917
sleep$date[which(as.Date(sleep$date)>as.Date("2018-12-01")) ] #values that were typo 
which(as.Date(sleep$date)>as.Date("2018-12-01")) #index of values that were typos 

# correcting the date values 
sleep$date[55:65]
sleep$date[59] = "2015-10-29"  # instead of "2105-10-29"

sleep$date[2220:2230]
sleep$date[2223] = "2017-09-14" # instead of "2917-09-14"

sleep$date[2705:2715]
sleep$date[2711] = "2017-11-15" # instead of "2917-11-15"

range(as.Date(sleep$date, format = "%Y-%m-%d"),na.rm=T) 

sum(is.na(sleep$date)) # date is missing for 145 obs 
sum(is.na(sleep$date[sleep$o2utv4h == 1])) # date is missing for ZERO obs in those with valid sat measurement 

# Subsetting the data for valid flow and sat monitoring ####
# Analysis using AHI should only include individuals with valid flow and sat monitoring 
# Analysis using ODI should only include individuals with valid sat monitoring 

#subsetting the dataset to only those with valid flow and sat monitoring 
valid.ahi = sleep[flutv4h==1 & o2utv4h==1,]

#subsetting the dataset to only those with valid sat monitoring (enough for analysis with ODI or Sat90%)
valid.odi = sleep[o2utv4h==1,]

# Summary Apnea-hipopnea index (AHI)
summary(valid.ahi$ahi)

# Summary Oxygen desaturation index (ODI)
summary(valid.odi$odi)

# Summary % of registration with sat≤90% 
summary(valid.odi$sat90)

# 2 individuals with missing data for AHI among those with valid flow and sat monitoring
valid.ahi[is.na(valid.ahi$ahi),c("idnr", "fldeutv_min", "spo2utv_min","date","anstart", "anstop", "ahi", "odi")]

# ODI missing for 4 individuals (# 2 are also missing informatio for AHI, 1 is also missing information on sat90%)
valid.odi$idnr[which(is.na(valid.odi$odi))]

# % of registration with sat≤90% missing for 1 individual 
valid.odi$idnr[which(is.na(valid.odi$sat90))]

# Creating a variable OSA severity 
valid.ahi$OSAcat = rec(valid.ahi$ahi, rec = "0:4.9=0; 5:14.9=1; 15:29.9=2 ; 30:max=3")
valid.ahi$OSAcat = factor(valid.ahi$OSAcat, levels = c(0,1,2,3), 
                          labels = c("no OSA", "Mild","Moderate", "Severe"))

#Frequency table by OSA cat 
tableOSAcat = freq(valid.ahi$OSAcat, style='rmarkdown')

#histogram of flow monitoring duration 
p1 = ggplot(data = sleep, aes(x=fldeutv_min)) + geom_histogram() + geom_vline(xintercept = 240, color = "red")  +
  scale_x_continuous(breaks = c(180, 240, 300, 420, 540, 660), labels = c("3h", "4h", "5h", "7h", "9h", "11h")) +
  ggtitle("Duration of flow monitoring") +  xlab("") +
  theme(axis.text.x = element_text(color = c("black", "red", "black", "black", "black", "black")))

#histogram of saturation monitoring duration 
p2 = ggplot(data = sleep, aes(x=spo2utv_min)) + geom_histogram() + geom_vline(xintercept = 240, color = "red")  +
  scale_x_continuous(breaks = c(180, 240, 300, 420, 540, 660), labels = c("3h", "4h", "5h", "7h", "9h", "11h")) +
  ggtitle("Duration of saturation monitoring") +  xlab("") +
  theme(axis.text.x = element_text(color = c("black", "red", "black", "black", "black", "black")))


# Histogram AHI
p3 = ggplot(data = valid.ahi, aes(x=ahi)) + geom_histogram()  +
  ggtitle("Hist Apnea-hypoapnea index") +  xlab("") + 
  geom_vline(xintercept = 5, color = "black", linetype = "twodash") + 
  geom_vline(xintercept = 15, color = "black", linetype = "twodash") + 
  geom_vline(xintercept = 30, color = "black", linetype = "twodash")

# Histogram ODI
p4 = ggplot(data = valid.odi, aes(x=odi)) + geom_histogram()  +
  ggtitle("Hist Oxygen desaturation index") +  xlab("")

# Histogram Sat90% 
p5 = ggplot(data = valid.odi, aes(x=sat90)) + 
  geom_histogram(color="black", fill="lightskyblue2") +
  geom_density(alpha=.2, fill="#FF6666") +
  ggtitle("Hist - T90%") +  xlab("") +
  theme_light()+
  theme(plot.title= element_text(hjust = 0.5, size = 16,face = "bold"))

ggsave("hist.flow", plot = p1, device = "png", 
path = "/home/baldanzi/Sleep_apnea/Descriptive/")

ggsave("hist.saturation", plot = p2, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

ggsave("hist.ahi.png",plot = p3, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

ggsave("hist.odi.png",plot = p4, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

ggsave("hist.sat90.png",plot = p5, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

# Uploading phenotype data ####
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
pheno=fread("SCAPIS-DATA-PETITION-170-20210315.csv", header = T, na.strings=c("", "NA"))
# Restricting dataset to Uppsala participants only 
pheno = pheno[Site=="Site5", ]


# Uploading MGS data ####
data.MGS = fread("/home/baldanzi/Datasets/MGS/MGS_relative_abundance_4839_upp_4980_malmo.tsv")
# 4839 Uppsala participants have gut microbiome data

# Adjusting ids so that all have the same format ####
ids=fread("/home/baldanzi/Datasets/scapis_idkey/id_conversion.txt", header=T) 
ids$id = as.character(ids$export_id)
nrow(ids[ids$subject_id %in% grep("5-",ids$subject_id, value = T ),] )
# 5033 Uppsala participants


data.MGS = merge(data.MGS, ids, by = "id", all=T, all.x=T, all.y=F)
setnames(data.MGS,"subject_id", "SCAPISid")

# Restricting to Uppsala participants 
data.MGS = data.MGS[data.MGS$SCAPISid %in% grep("5-",data.MGS$SCAPISid, value = T ),]

# Place of birth data ####
placebirth = fread("/home/baldanzi/Datasets/Placeofbirth.csv", header=T, sep=",")
placebirth$pob = factor(placebirth$q005a, 
                        levels = c("scandinavia", "europe", "asia", "other"), 
                        labels = c("Scandinavia", "Europe", "Asia", "other"))
placebirth[,q005a:=NULL]
#restricting to Uppsala 
placebirth = placebirth[placebirth$SCAPISid %in% grep("5-",placebirth$SCAPISid, value = T ),]

# Diet data #### (the data information in use is also available at the _pheno_ file)
#diet = fread("/home/baldanzi/Datasets/diet_SCAPIS/export-petition-122-V0.0.50-2021-01-13.csv", header=T)
#a = c("Subject", "Fibrer")
#diet = diet[,a,with=F]
#setnames(diet, "Subject", "SCAPISid")

# Metabolon data ####
metabolon=fread("/home/baldanzi/Datasets/Metabolon/clean/scapis_metabolon_scapisid.tsv", header=T)

# Metabolon Metadata - collection date ####
# Import metabolon metadata 

metab_metadata=fread('/home/baldanzi/Datasets/Metabolon/clean/metabolon_metadata_scapisid.tsv', header=T)
a = c("SCAPISid", "COLLECTION_DATE")
metab_collection_date = metab_metadata[,a,with=F]
setnames(metab_collection_date,"COLLECTION_DATE", "metabolon_collection_date")

# 28 participants have sleep and gut microbiota data but do not have metabolomics data


# Recoding variables ####

#Age
pheno[,age:=AgeAtVisitOne]

#Sex
pheno$Sex = factor(pheno$Sex, levels = c("MALE", "FEMALE"), 
                       labels = c("male", "female"))

#smoking 
setnames(pheno, 'derived_smoke_status', 'smokestatus')
pheno$smokestatus = rec(pheno$smokestatus, 
                        rec = "NEVER=never; CURRENT=current; EX_SMOKER=former; NA,UNKNOWN=NA")
pheno$smokestatus = factor(pheno$smokestatus, levels = c("never", "former", "current"))

# Highest education 
pheno$educat = rec(pheno$cqed001, rec = "-99=NA; else=copy")
pheno$educat = factor(pheno$educat, levels = c(0,1,2,3), 
                      labels = c("uncompleted primary or lower secondary", 
                                 "lower secondary education", "upper secondary education", 
                                 "university education"))

# BMI categories 
pheno$BMIcat = rec(pheno$BMI, rec = "min:24.9=1; 25:29.9=2; 30:max = 3")
pheno$BMIcat = factor(pheno$BMIcat, levels = c(1,2,3), labels = c("<25","25-30",">=30"))

# Self-reported physical activity 
pheno$leisurePA = rec(pheno$cqpa012, rec = "4=NA; else=copy")
pheno$leisurePA = factor(pheno$leisurePA, levels = c(0,1,2,3), 
                        labels = c("mostly sedentary", "moderate activity", 
                                   "regular and moderate activity ", "regular exercise or training"))

# Diabetes variable
pheno$diabd = rec(pheno$Diabetes, rec = "NORMOGLYCEMIA=0; ELEV_HBA1C,IFG=1; KNOWN_DM,NEW_DM=2")
pheno$diabd = factor(pheno$diabd, levels = c(0,1,2), 
                       labels = c("normoglycemic", "impaired glucose tolerance", "type 2 diabetes" ))

# Fiber adjusted for energy intake ####
# Removing over- and under- reported based on the 3SD of the natural log of Energy intake
pheno[Energi_kcal<200 | Energi_kcal >8000, Energi_kcal:=NA] # Removing invalid values
pheno[,log.energi:= log(pheno$Energi_kcal)] # Calculating the log energy intake

# Estimate sd and mean of energy intake
male.sd <-  sd(pheno[Sex=="male",log.energi],na.rm=T)
male.mean <- mean(pheno[Sex=="male",log.energi],na.rm=T)
female.sd <- sd(pheno[Sex=="female",log.energi],na.rm=T)
female.mean <- mean(pheno[Sex=="female",log.energi],na.rm=T)

# Using 3 sd to categorize over- and under-reporters
pheno[!is.na(Energi_kcal),energi.reporter:="ok"]

ll <- female.mean-(3*female.sd)
ul <- female.mean+(3*female.sd)
pheno[Sex=="female" & log.energi<ll, energi.reporter:="under"]
pheno[Sex=="female" & log.energi>ul, energi.reporter:="over"]

ll <- male.mean-(3*male.sd)
ul <- male.mean+(3*male.sd)
pheno[Sex=="male" & log.energi<ll, energi.reporter:="under"]
pheno[Sex=="male" & log.energi>ul, energi.reporter:="over"]

# Removing under- or over- reporter
pheno[,energi.original:=Energi_kcal]
pheno[energi.reporter!="ok",Energi_kcal:=NA]
# Consider revising this variable  after excluding those with unreliable FFQ reporting. 
pheno[,fiber.kcal:=(Fibrer/Energi_kcal)*1000]


# Self-reported hypertension variable
pheno$hypertension = factor(pheno$cqhe034, levels = c("NO", "YES"), labels = c("no", "yes"))

# Self-reported dyslipidemia 
pheno$dyslipidemia = factor(pheno$cqhe036, levels = c("NO", "YES"), labels = c("no", "yes"))

# Self-reported diabetes medication 
pheno$diabmed = factor(pheno$cqhe039, levels = c("NO", "YES"), labels = c("no", "yes"))

# Self-reported hypertension medication 
pheno$hypermed = factor(pheno$cqhe035,  levels = c("NO", "YES"), labels = c("no", "yes"))

# Self-reported dyslipidemia medication 
pheno$dyslipmed = factor(pheno$cqhe037,  levels = c("NO", "YES"), labels = c("no", "yes"))

# Self-reported Sleep apnea diagnosis 
pheno$apnea_self = factor(pheno$cqhe058,  levels = c("NO", "YES"), labels = c("no", "yes"))

# Self-reported Sleep apnea treatment  
pheno$apneatto_self = factor(pheno$cqhe060,  levels = c("NO", "YES"), labels = c("no", "yes"))

# Self-reported Sleep apnea treatment CPAP
pheno$cpap_self = ifelse(pheno$cqhe061=="CPAP", "yes", "no")

# Self-reported Sleep apnea treatment Oral appliance 
pheno$splint_self = ifelse(pheno$cqhe061=="SPLINT", "yes", "no")

# Self-reported Sleep apnea treatment Surgery 
pheno$apneasurgery_self = ifelse(pheno$cqhe061=="SURGERY", "yes", "no")

# CPAP and Splint data 
# in the original data, CPAP and Splint only have the values 1 and NA
valid.ahi$cpap = factor(valid.ahi$cpap, levels = c(0,1), labels = c("no", "yes"))
valid.ahi$splint = factor(valid.ahi$splint, levels = c(0,1), labels = c("no", "yes"))
valid.ahi$cpap[is.na(valid.ahi$cpap)] = "no"
valid.ahi$splint[is.na(valid.ahi$splint)] = "no"

valid.odi$cpap = factor(valid.odi$cpap, levels = c(0,1), labels = c("no", "yes"))
valid.odi$splint = factor(valid.odi$splint, levels = c(0,1), labels = c("no", "yes"))
valid.odi$cpap[is.na(valid.odi$cpap)] = "no"
valid.odi$splint[is.na(valid.odi$splint)] = "no"

# Medication data based on metabolomics measurement
a = c("SCAPISid","MET_100002725","MET_100002808","MET_100002405")
medication = metabolon[,a,with=F]
medication[,ppi:=ifelse(MET_100002725>min(MET_100002725,na.rm=T) |
                      MET_100002808>min(MET_100002808,na.rm=T),1,0)]
medication[,metformin:=ifelse(MET_100002405>min(MET_100002405,na.rm=T),1, 0)]

medication[is.na(ppi), ppi:=0]
medication[,ppi:=factor(ppi, levels = c(0,1), 
                 labels = c("no", "yes"))]
medication[is.na(metformin), metformin:=0]
medication[,metformin:=factor(metformin, levels = c(0,1), 
                        labels = c("no", "yes"))]

a = c("SCAPISid","ppi","metformin")
medication = medication[,a,with=F]

# Epworth Sleepiness Scale score ####
# varibles for the ESS score
essvar = c("cqsh006", "cqsh007","cqsh008","cqsh009","cqsh010","cqsh011","cqsh012","cqsh013")
#Transforming "do not want to answer" to NA
temp = apply(pheno[,essvar,with=F], 2, function(x){
  ifelse(x==-99,NA,x)
})
pheno[,(essvar):=as.data.table(temp)] 
# Summing up the ESS variables to create the ESS score
pheno[,ESS:=as.data.table(apply(pheno[,essvar,with=F], 1,sum))]


# Merging pheno (phenotype data) with valid.ahi data (data valid for analysis with AHI) ####
nrow(valid.ahi) # 3301 individuals with valid flow and sat monitoring 
nrow(pheno) # 5036 Uppsala individuals with phenotype data
setnames(valid.ahi, "id", "SCAPISid") #renaming the id variable
setnames(valid.odi, "id", "SCAPISid") #renaming the id variable
setnames(pheno, "Subject", "SCAPISid") #renaming the id variable

#vali.ahi + pheno
sum(valid.ahi$SCAPISid %in% pheno$SCAPISid) # 3301 with valid monitoring and pheno data 
valid.ahi = merge(valid.ahi, pheno, by = "SCAPISid", all.x=T, all.y=F )

#valid.ahi + pheno + data.MGS
sum(valid.ahi$SCAPISid %in% data.MGS$SCAPISid) # 3208 with valid monitoring and MGS data 
valid.ahi = merge(valid.ahi, data.MGS, by = "SCAPISid", all=T, all.x=F, all.y=F)

#valid.ahi + pheno + data.MGS + place of birth
valid.ahi = merge(valid.ahi, placebirth, by = "SCAPISid", all.x=T, all.y=F)

#valid.ahi + pheno + data.MGS + place of birth + diet 
#valid.ahi = merge(valid.ahi, diet, by ="SCAPISid", all.x=T, all.y=F)

# ... + medication 
valid.ahi = merge(valid.ahi,medication, by="SCAPISid", all.x=T, all.y=F)

# ... + metabolon collection date
valid.ahi=merge(valid.ahi, metab_collection_date, by="SCAPISid", all.x=T, all.y=F)

# 28 participants have sleep and gut microbiota data but do not have metabolomics data
sum(valid.ahi[,SCAPISid] %in% metabolon[,SCAPISid])
nrow(valid.ahi)-sum(valid.ahi[,SCAPISid] %in% metabolon[,SCAPISid])

# Changing order
setcolorder(valid.ahi, c("SCAPISid", "age", "Sex"))

# Excluding 2 individuals with missing information on AHI
valid.ahi = valid.ahi[!is.na(ahi),]

# Changing order
setcolorder(valid.ahi, c("SCAPISid", "age", "Sex"))

# Excluding 2 individuals with missing information on AHI
valid.ahi = valid.ahi[!is.na(ahi),]


# Clinical microbiomics variables ####
cmvar = fread("/home/baldanzi/Datasets/MGS/CM_lab_variables_4838_upp_4980_malmo.tsv",header=T, na.strings=c("", "NA")) 
selected_id = fread("/home/baldanzi/Datasets/scapis_idkey/selected_participants_id_sample_id.tsv", 
                    header=T, na.strings=c("", "NA"))
cmvar = merge(selected_id, cmvar, by = "sample.id", all=F, all.x=T )
ids=fread("/home/baldanzi/Datasets/scapis_idkey/id_conversion.txt", header=T) 
ids$id = as.character(ids$export_id)
cmvar = merge(cmvar, ids, by = "id", all=F, all.x=T)
cmvar[is.na(subject_id), subject_id:=id]
setnames(cmvar, "subject_id", "SCAPISid")
setnames(cmvar,"dna.yield..\302\265g.", "dna.yield")

#merging CM variables with the valid AHI data. 
a = c("SCAPISid", "plate", "received", "read.pairs.per.sample", "dna.yield")
valid.ahi = merge(valid.ahi, cmvar[,a,with=F], by="SCAPISid", all=F, all.x=T, all.y=F)
valid.ahi$plate = as.factor(valid.ahi$plate)
valid.ahi$received = as.factor(valid.ahi$received)

#Saving the data #### 
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
fwrite(valid.ahi, file = "validsleep.MGS.Upp.tsv", sep = "\t")

# Saving those variables that were managed in Rds format 
dat1 = valid.ahi %>% select(SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, WaistHip, educat,
                            leisurePA, pob, diabd, hypertension, dyslipidemia, diabmed, 
                            hypermed, dyslipmed, ppi, Fibrer, Energi_kcal,ESS,apnea_self, apneatto_self,
                            cpap_self, splint_self, apneasurgery_self,
                            ahi, odi, sat90, cpap, splint)
save(dat1, file = "data_table1")



# Nicotine metabolites ("cotinine", "hydroxycotinine", "cotinine N-oxide") ####
a = c("SCAPISid" ,"MET_848","MET_100002717", "MET_100002719")
nicotine=metabolon[,a,with=F]

# Merging nicotine metabolites with smokestatus 
a = c("SCAPISid", "smokestatus", "cqto015")
nicotine = merge(nicotine, valid.ahi[,a, with=F], by="SCAPISid", all.x = F, all.y=T)
setnames(nicotine,"cqto015", "snus_ever1mo")

# Saving nicotine metabolites data 
setwd("/home/baldanzi/Datasets/Metabolon")
fwrite(nicotine, file = "nicotine_metab.csv", sep=",")

# Saving valid.odi data 
valid.odi = merge(valid.odi, pheno, by = "SCAPISid", all.x=T, all.y=F )
valid.odi = merge(valid.odi, data.MGS, by = "SCAPISid", all=T, all.x=F, all.y=F)
valid.odi = merge(valid.odi, placebirth, by = "SCAPISid", all.x=T, all.y=F)
a = c("SCAPISid", "plate", "received", "read.pairs.per.sample", "dna.yield")
valid.odi = merge(valid.odi, cmvar[,a,with=F], by="SCAPISid", all=F, all.x=T, all.y=F)
valid.odi = merge(valid.odi,medication, by="SCAPISid", all.x=T, all.y=F)
fwrite(valid.odi, file = "/home/baldanzi/Datasets/sleep_SCAPIS/validodi.MGS.Upp.tsv", sep = "\t")

# 3622 individuals have a valid T90 measurement 

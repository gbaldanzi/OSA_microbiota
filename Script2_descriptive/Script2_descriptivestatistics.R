# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# Last update: 2021-09-22

message("#Descriptive Statistics#")
print(Sys.time())

# Loading packages
library(tidyverse)
library(data.table)
library(summarytools)
library(ggpubr)
library(Hmisc)
library(compareGroups)
library(flextable)
library(vegan)
library(RColorBrewer)


  # Import data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

#Variables are: SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, WaistHip, educat,
#leisurePA, placebirth, diabd, hypertension, dyslipidemia, diabmed, 
#hypermed, dyslipmed, ppi, fiber, EI,ESS,apnea_self, apneatto_self,
#cpap_self, splint_self, apneasurgery_self,
#ahi, odi, sat90, cpap, splint 
  
  dat1 = pheno[valid.ahi=='yes',] %>% select(SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, WaistHip, educat,
                                             leisurePA, placebirth, diabd, hypertension, dyslipidemia, diabmed, 
                                             hypermed, dyslipmed, ppi, Fibrer, Energi_kcal,ESS,apnea_self, apneatto_self,
                                             cpap_self, splint_self, apneasurgery_self,
                                             ahi, odi, sat90, cpap, splint, shannon)

# Variable diagnostics ####
# Number of individuals 
nrow(dat1) #3206


# Checking for missingness in the variables. 
miss = data.frame(Missing = apply(dat1,2,function(x){sum(is.na(x))}))
a = apply(dat1,2,function(x){sum(is.na(x))/length(x)})
miss = cbind(Variables = rownames(miss), miss, Percentage = paste0(round(a*100,0)," %"))
#flextable(miss)

# Number of individuals having at least one missing value for the covariates
#perrow = apply(dat1, 1, function(x) {any(is.na(x))})
#sum(perrow) 

# Summary of variables
summary = sapply(names(dat1), function(x) {summary(dat1[[x]])})

# "Table 1" - Population characteristics ####
datatable1 = as.data.frame(dat1[,-1])  # Passing data to a new object without the SCAPIS ID

# For the purpose of this table, individuals with missing values for the following variables
# were included in the "no" category. (only category "yes" is displayed in the table)
for(i in which(names(datatable1) %in% c('hypertension', 'dyslipidemia', 'diabmed',
                                'hypermed','dyslipmed', 'ppi','apnea_self', 'apneatto_self',
                                'cpap_self', 'splint_self', 'apneasurgery_self','cpap', 'splint'))){
 datatable1[is.na(datatable1[i]),i] = "no"
}

#Creating labels for the variables 
mylabel <- c("OSA severity", "Age (yrs)", "Sex", "Smoking status", "Alcohol intake (g)", 
             "BMI (kg/m2)", "WHR","Highest education achieved", "Self-reported leisure physical activity", 
             "Place of birth", "Type 2 diabetes", "Hypertension", "Dyslipidemia",
             "Diabetes medication",
             "Hypertension medication", "Dyslipidemia medication", "PPI","Fiber (g/day)","Energy intake",
             "ESS", "Sleep apnea (self-reported)", "Sleep apnea treatment", "CPAP", "Oral appliance",
             "Surgery", "AHI (events/h)", "ODI (events/h)", 
             "T90%", "CPAP during examination", "Oral appliance during examination", "Shannon Index")

j=1
for(i in names(datatable1)){
  label(datatable1[[i]]) <- mylabel[j]
  j=j+1
}

# Create table of population characteristics by OSA group (no OSA, Mild, Moderate, or Severe OSA)
t = compareGroups(OSAcat ~ age + Sex + smokestatus + Alkohol + BMI + educat +
                    leisurePA + placebirth + diabd + hypertension + dyslipidemia + diabmed + 
                    hypermed + dyslipmed+ ppi+Fibrer+Energi_kcal+ESS+apnea_self+apneatto_self+
                    cpap_self+ splint_self+ apneasurgery_self+ahi+ odi+ sat90+
                    ESS+ cpap+ splint + shannon, data= datatable1, 
                  include.miss = FALSE, chisq.test.perm = TRUE)
t1 = createTable(t, hide.no = "no")
t1

saveRDS(t1, file='/home/baldanzi/Sleep_apnea/Descriptive/sleepapnea_table1.rds')

rm(datatable1)

#______________________________________________________________________________

# Descriptive statistics ####

  valid.ahi <- pheno[valid.ahi=='yes',]

#---------------------------------------------------------------------------#
# We wanted to check if blood samples for metabolomics and Sleep records (ApneaLink)
#took place close to each other (difference in calender data between the two)


#### Dates ####
# Was sleep apnea assessed close to when metabolomics was collected? 
valid.ahi[,metabolon_collection_date := as.Date(metabolon_collection_date, 
                                              format = "%d-%b-%y")]
a = c("date", "metabolon_collection_date", "AnthropometryCollectionDate","LdlSamplingDate" )
# Number of participants for which metabolon_collection_date is missing 
nrow(valid.ahi[is.na(metabolon_collection_date),]) # 1229

# Number of participants for which metabolon_colletion and sleep apnea assessment date 
# are the same
nrow(valid.ahi[date==metabolon_collection_date,]) # 1502

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
            "diabd", "placebirth", "Alkohol", "Fibrer", "Energi_kcal")

# Selecting 
contain_missing = which(apply(dat1[,listvar,with=F], 1, function(x){any(is.na(x))}))
dat1[,incomplete_obs:=NULL]
dat1[contain_missing,incomplete_obs:=1]
dat1[!contain_missing,incomplete_obs:=2]
dat1$incomplete_obs = factor(dat1[,incomplete_obs], 
                               levels= c(1,2),
                               labels = c("Incomplete obs.", 
                                          "Complete obs."))


# Preparing data for table: Participants characteristics divided by 
#observation complete or incomplete 
  datafortable = as.data.frame(dat1[,-1]) 

# Labelling variables 
  mylabel <- c("OSA severity", "Age (yrs)", "Sex", "Smoking status", "Alcohol intake (g)", 
             "BMI (kg/m2)", "WHR","Highest education achieved", "Self-reported leisure physical activity", 
             "Place of birth", "Type 2 diabetes", "Hypertension", "Dyslipidemia",
             "Diabetes medication",
             "Hypertension medication", "Dyslipidemia medication", "PPI","Fiber (g/day)",
             "Energy intake","ESS",
             "Sleep apnea (self-reported)", "Sleep apnea treatment", "CPAP", "Oral appliance",
             "Surgery",
             "AHI (events/h)", "ODI (events/h)", 
             "T90%", "CPAP during examination", "Oral appliance during examination","Shannon diversity",
             "Incomplete observation" )
  
  j=1
  for(i in names(datafortable)){
     label(datafortable[[i]]) <- mylabel[j]
     j=j+1
  }

# Creating the table 
t = compareGroups(incomplete_obs ~ ahi + shannon + age + Sex + BMI+  Alkohol + Fibrer + Energi_kcal + smokestatus + educat +
                    leisurePA + placebirth + diabd + hypertension+apnea_self, data= datafortable, 
                  include.miss = FALSE, chisq.test.perm = TRUE)
t1 = createTable(t)
t1

saveRDS(t1, file='/home/baldanzi/Sleep_apnea/Descriptive/completVSincompletObs.rds')

#---------------------------------------------------------------------------#
#### Clinical microbiomics variables ####
# Variables related to how fecal samples were analysed to produce the gut microbiota data

# 
#tcm = compareGroups(OSAcat~read.pairs.per.sample + dna.yield, data= valid.ahi)
#t1 = createTable(tcm)

#---------------------------------------------------------------------------#

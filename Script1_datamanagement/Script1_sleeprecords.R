# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2021-09-30

# Sleep records were an add-on to the SCAPIS study in Uppsala when participants
# were assessed for sleep apnea during one night's sleep with the ApneaLink air flow 
# and oxygen saturation devide. These data were provided to us by Eva Lindberg at from 
# the Respiratory-, allergy- and sleep research group

# This script will perform minor adjustments to the sleep data. Importantly, it
# creates two flag variables: valid.ahi and valid.t90

## valid.ahi indicates participants to be included in analysis with the AHI variable
## valid.t90 indicates participants to be included in analysis with the variables T90 and ODI

# Data management ####

# Loading packages
pacman::p_load(tidyverse, chron, Hmisc, sjmisc, data.table)

rm(list=ls())

# Import sleep recording data 
sleep <- fread('/home/baldanzi/Datasets/sleep_SCAPIS/original_data/Sleep_SCAPIS_Uppsala_Final.csv',
               na.strings = c("","NA",NA))
sleep <- sleep[AttendSleepStudy==TRUE,]

setnames(sleep, "id", "SCAPISid") #renaming the id variable

# Data summary 
attach(sleep)
nr.idnr.valid = nrow(sleep[BothFlO2utv4h==1,]) # number of ind with at least 4h of flow and sat monitoring 
fl = summary(fldeutv_min) 
fl = substr(times((fl%/%60 +  fl%%60 /60)/24), 1, 5) # Mean and median duration of flow monitoring 

sat = summary(spo2utv_min)
sat = substr(times((sat%/%60 +  sat%%60 /60)/24), 1, 5) # Mean and median duration of flow monitoring
detach(sleep)

# Date of recording ####

# Date interval 
sleep[,Date := as.Date(date,format= c("%d%b%Y"))]
sleep[is.na(Date),Date := as.Date(sleep$date[is.na(sleep$Date)],format= c("%d-%b-%y"))]

range(sleep$Date,na.rm=T) #highest value with year = 2105
sleep$Date[which(sleep$Date>as.Date("2018-12-01")) ] #typo: "2105-10-29" 
which(sleep$Date>as.Date("2018-12-01")) #index of values that were typos 
  
  # Correct typo 
  sleep$Date[55:65]
  sleep$Date[59] <-  "2015-10-29"  # instead of "2105-10-29"
  
  range(sleep$Date,na.rm=T) 
  
  sum(is.na(sleep$Date)) # date is missing for 145 obs 
  sum(is.na(sleep$Date[sleep$o2utv4h == 1])) # date is missing for ZERO obs in those with valid sat measurement 
  

# Subsetting the data for valid flow and sat monitoring ####
  
  # Analysis using AHI should only include individuals with valid flow and sat monitoring 
  # Analysis using T90 should only include individuals with valid sat monitoring 

#subsetting the dataset to only those with valid flow and sat monitoring 
  sleep[,valid.ahi:=ifelse(sleep$BothFlO2utv4h==1, "yes", "no")]
  sleep[,valid.t90:=ifelse(sleep$o2utv4h==1, "yes", "no")]



# 2 individuals with missing data for AHI among those with valid flow and sat monitoring
  sleep[is.na(ahi) & valid.ahi=="yes",c("SCAPISid", "fldeutv_min","Date", "ahi", "sat90")]

# T90 missing for 1 individual 
  sleep[is.na(sleep$sat90) & valid.t90=='yes',.(SCAPISid,sat90)]

# Excluding 2 individuals with missing information on AHI
  sleep[valid.ahi=='yes' & is.na(ahi),valid.ahi:='no']

# Excluding 2 individuals with missing information on T90
  sleep[valid.t90=='yes' & is.na(sat90),valid.t90:='no']

# Creating a variable OSA severity 
sleep[valid.ahi=="yes",OSAcat:= rec(ahi, rec = "0:4.9=0; 5:14.9=1; 15:29.9=2 ; 30:max=3")]
sleep[,OSAcat:= factor(OSAcat, levels = c(0,1,2,3), 
                          labels = c("no OSA", "Mild","Moderate", "Severe"))]

# CPAP and Splint data 
# in the original data, CPAP and Splint only have the values 1 and NA
sleep[,cpap:= factor(sleep$cpap, levels = c(0,1), labels = c("no", "yes"))]
sleep[,splint:= factor(sleep$splint, levels = c(0,1), labels = c("no", "yes"))]
sleep[is.na(sleep$cpap),cpap:= "no"]
sleep[is.na(sleep$splint),splint:= "no"]
sleep[,t90:= sat90]


saveRDS(sleep,"/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording/sleep.rds")


  
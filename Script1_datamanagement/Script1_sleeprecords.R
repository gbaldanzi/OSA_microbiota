# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script for data management of sleeping records

# Data management ####

# Loading packages
pacman::p_load(tidyverse, grid, chron, rio, Hmisc, sjmisc, summarytools, data.table)

rm(list=ls())

# Import sleep recording data 

setwd("/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording")
sleep = fread("Data_Sleep_SCAPIS_Uppsala.csv", header=T, na.strings=c("", "NA"))
setnames(sleep, "id", "SCAPISid") #renaming the id variable


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

  # sleep duration during the examination 

  to.min.fun = function(a){
    
    splitted = strsplit(as.character(a),':')
    in.min = sapply(splitted,function(x){
      as.numeric(x[1])*60 + as.numeric(x[2])
    })
    return(in.min)
  }
    
    
  sleep$sovtid <- sub("8.","8:",sleep$sovtid)
  sleep[,sovtid:=times(paste0(sovtid,':00'))]
  
  sleep[,sleeptime:=to.min.fun(sovtid)]

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
valid.t90 = sleep[o2utv4h==1,]

# 2 individuals with missing data for AHI among those with valid flow and sat monitoring
valid.ahi[is.na(valid.ahi$ahi),c("idnr", "fldeutv_min", "spo2utv_min","date","anstart", "anstop", "ahi", "odi")]

# ODI missing for 4 individuals (# 2 are also missing informatio for AHI, 1 is also missing information on sat90%)
valid.t90$idnr[which(is.na(valid.t90$odi))]

# % of registration with sat≤90% missing for 1 individual 
valid.t90$idnr[which(is.na(valid.t90$sat90))]

# Creating a variable OSA severity 
valid.ahi$OSAcat = rec(valid.ahi$ahi, rec = "0:4.9=0; 5:14.9=1; 15:29.9=2 ; 30:max=3")
valid.ahi$OSAcat = factor(valid.ahi$OSAcat, levels = c(0,1,2,3), 
                          labels = c("no OSA", "Mild","Moderate", "Severe"))

# CPAP and Splint data 
# in the original data, CPAP and Splint only have the values 1 and NA
valid.ahi$cpap = factor(valid.ahi$cpap, levels = c(0,1), labels = c("no", "yes"))
valid.ahi$splint = factor(valid.ahi$splint, levels = c(0,1), labels = c("no", "yes"))
valid.ahi$cpap[is.na(valid.ahi$cpap)] = "no"
valid.ahi$splint[is.na(valid.ahi$splint)] = "no"

valid.t90$cpap = factor(valid.t90$cpap, levels = c(0,1), labels = c("no", "yes"))
valid.t90$splint = factor(valid.t90$splint, levels = c(0,1), labels = c("no", "yes"))
valid.t90$cpap[is.na(valid.t90$cpap)] = "no"
valid.t90$splint[is.na(valid.t90$splint)] = "no"

valid.t90$t90 = valid.t90$sat90


  # Saving the datasets ####
  saveRDS(valid.ahi,"/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording/valid.ahi.rds")
  saveRDS(valid.t90,"/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording/valid.t90.rds")
  

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
p4 = ggplot(data = valid.t90, aes(x=odi)) + geom_histogram()  +
  ggtitle("Hist Oxygen desaturation index") +  xlab("")

# Histogram Sat90% 
p5 = ggplot(data = valid.t90, aes(x=sat90)) + 
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

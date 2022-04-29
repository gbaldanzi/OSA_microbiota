# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script for preparing the data. 

# This script will prepare the data for the analysis on the association between 
# obstructive sleep apnea and gut microbiota features. 

# Prepare here means creating new variables and exclude the diet variables for under- and 
# over-reporters of diet intake. It also removes species whose prevalence is <1% and 
# calculates the Shannon diversity index. 

  # Packages
  library(tidyverse)
  library(Hmisc)
  library(sjmisc)
  library(data.table) 
  library(vegan)

  # input = folder containing the data 
  input <- "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/data_processed/"
  
  # work folder to save the final data
  #if(!dir.exists("work")){dir.create("work")}
  # work <- './work/'
  work <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'

  # results folder 
  #if(!dir.exists("results")){dir.create("results")}
  # results <- './results/'
  
  # Import data 
  pheno <-readRDS(paste0(input,"pheno_sleep_mgs.rds"))


  # Prepare variables ####

  # Month of anthropometric collection date - variable to account for season
  pheno[,visit.month:=format(as.POSIXct(pheno$AnthropometryCollectionDate),"%B")]
  pheno[,visit.month:=factor(visit.month,c(month.name,"June.July"))]
  
      # Merging June to July due to the low number of July participants 
      pheno[visit.month %in% c("June","July"),visit.month:="June.July"]
      
  
  # Excluding diet information for over- and under- reporters based on the 3SD of the natural log of Energy intake

  # remove extreme outlier
  pheno[Energi_kcal==0, Energi_kcal:=NA]
  
  # Calculating the log energy intake
  pheno[,log.energi:= log(pheno$Energi_kcal)] 

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

  
  # Removing diet data from under- or over- reporter
  pheno[,energi.original:=Energi_kcal]
  pheno[energi.reporter!="ok",Energi_kcal:=NA]
  pheno[energi.reporter!="ok",Fibrer:=NA]
   
 
  #Those with valid AHI values have at least 4 hours of air flow and oxygen saturation monitoring
  pheno[!is.na(ahi), valid.ahi := ifelse(BothFlO2utv4h==1, "yes", "no")]
  
  #Those with valid T90 and ODI values have at least 4 hours of oxygen saturation monitoring
  pheno[!is.na(t90) & !is.na(odi), valid.t90 := ifelse(o2utv4h==1, "yes", "no")]
  

  # Exclude CPAP users 
  pheno[cqhe061 == 'CPAP', c("valid.ahi","valid.t90") := "no"] # Self-reported CPAP use for OSA (questionnaire)
  pheno[cpap == 'yes', c("valid.ahi","valid.t90") := "no"] # Self-reported CPAP use during OSA assessment
  
  
  # Assigning to NA the AHI, T90 and ODI values for CPAP users or  with monitoring < 4 hours
  pheno[valid.ahi=='no', ahi:=NA]
  pheno[valid.t90=='no', c("t90","odi"):=NA]
  
  
  # Create groups of OSA severity based on AHI
  pheno[valid.ahi=="yes", OSAcat:= rec(ahi, rec = "0:4.9=0; 5:14.9=1; 15:29.9=2 ; 30:max=3")]
  pheno[,OSAcat:= factor(OSAcat, levels = c(0,1,2,3), 
                         labels = c("No OSA", "Mild","Moderate", "Severe"))]
  
  
  # Create groups of OSA severity based on T90
  pheno[t90==0, t90cat := 'T90 = 0' ]
  pheno[t90!=0, t90cat :=  cut(t90,breaks = quantile(t90,
                                                     probs = seq(0,1,by=1/3),
                                                     na.rm=T), include.lowest = T)]
 
  pheno[,t90cat:=factor(t90cat, levels = c("T90 = 0", "[1,3]", "(3,14]","(14,100]" ),
                        labels = c("T90 = 0","t1","t2","t3"))]
  
  
  # Create groups of OSA severity based on ODI
  pheno[,odicat := as.factor( cut(pheno$odi,breaks = quantile(odi, probs = seq(0,1,by=.25), na.rm=T), 
                                  include.lowest = T) )]
  pheno[,odicat := factor(odicat, levels = levels(odicat), labels = c("q1", "q2", "q3", "q4"))]

  
  
  # Calculate Shannon diversity index ####
  
  # Removing the MGS that are present in less than 1% of the individuals (prevalence<1%). 
  
  # Calculating MGS prevalence 
  species.names_all <- grep("HG3A",names(pheno),value=T)
  # presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
  data_pa <- decostand(x = pheno[,species.names_all, with=F], "pa")
  # calculate sum per species
  data_sum <- data.frame(prevalence = apply(data_pa, 2, sum)/nrow(data_pa))
  data_sum$MGS = rownames(data_sum)
  
  low.prevalence.species = data_sum$MGS[data_sum$prevalence <= 1/100] # 1,602 species have a prevalence greater than 1%
  pheno <- pheno[ , -low.prevalence.species, with=F] 
  
  # Calculate Shannon diversity 
  species.names_1perc = grep("HG3A",names(pheno),value=T) # vector with MGS names 
  pheno[, shannon := diversity(pheno[,species.names_1perc, with=F],index="shannon")]
  
  
  # Save data
  saveRDS(pheno, file=paste0(work,"pheno_sleep_mgs_shannon.rds"))
  

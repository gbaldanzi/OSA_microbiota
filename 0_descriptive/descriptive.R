# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# Last update: 2022-02-01

# Script to create a Table 1 with participants baseline characteristics 

# Loading packages
library(tidyverse)
library(data.table)
library(Hmisc)
library(compareGroups)

  # Folders 
  descriptive.folder <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/descriptive/'
  work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'

  # Import data
  pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))


# Table 1 only includes those participants with valid AHI 

    #Variables are: SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, WaistHip, educat,
  #leisurePA, placebirth, diabd, hypertension, dyslipidemia, diabmed, 
  #hypermed, dyslipmed, ppi, fiber, EI,ESS,apnea_self, apneatto_self,
  #cpap_self, splint_self, apneasurgery_self,
  #ahi, odi, sat90, cpap, splint 
  
  dat1 = pheno %>% 
    select(SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, educat,
                                           leisurePA, placebirth, diabd, hypertension, dyslipidemia, metformin, 
                                           hypermed, dyslipmed, ppi, Fibrer, Energi_kcal, ahi, odi,t90, shannon,
                                          t90cat, odicat, valid.t90, valid.ahi)


  dat1[OSAcat == "no OSA", OSAcat:= "No OSA"]
  dat1[,OSAcat:=factor(OSAcat, levels = c("No OSA", "Mild", "Moderate", "Severe"))]
  
  # "Table 1" - Population characteristics ####
  datatable1 <-  as.data.frame(dat1)  # Passing data to a new object without the SCAPIS ID

  # For the purpose of this table, individuals with missing values for the following variables
  # were included in the "no" category. (only category "yes" is displayed in the table)
  
  for(i in which(names(datatable1) %in% c('hypertension', 'dyslipidemia', 'metformin',
                                        'hypermed','dyslipmed', 'ppi'))){
    datatable1[is.na(datatable1[i]),i] = "no"
  }

  
  #Creating labels for the variables 
  
  mylabel <- c("SCAPISid","OSA severity", "Age (yrs)", "Sex", "Smoking status", "Alcohol intake (g/day)", 
             "BMI (kg/m2)", "Education", "Leisure physical activity", 
             "Birth place", "Type 2 diabetes", "Hypertension", "Dyslipidemia",
             "Metformin",
             "Anti-hypertensive med.", "Hyperlipidemia med.", "PPI","Fiber (g/day)","Energy intake (kcal/day)",
             "AHI (events/h)","ODI (events/h)", "T90 (%)", "Shannon Index", "T90 groups", "ODI groups")

  j=1
  for(i in names(mylabel)){
    label(datatable1[[i]]) <- mylabel[j]
    j=j+1
  }

  # Create table of population characteristics by OSA group (no OSA, Mild, Moderate, or Severe OSA)
  t = compareGroups(OSAcat ~ age + Sex + smokestatus + Alkohol + BMI + educat +
                    leisurePA + placebirth + diabd + hypertension + dyslipidemia + metformin + 
                    hypermed + dyslipmed+ ppi+Fibrer+Energi_kcal+ahi+ odi+ t90+ shannon, 
                  data= datatable1, 
                  include.miss = FALSE, chisq.test.perm = TRUE, 
                  method = c(age = 2, Alkohol = 2, BMI = 2, 
                             ahi = 2 , odi = 2, t90 = 2, 
                             shannon = 2, Fibrer = 2, Energi_kcal = 2))
  
  t1 <-  createTable(t, hide.no = "no", show.p.overall = FALSE, show.all=T, 
                     digits = c(age = 1, Sex = 1, Alkohol = 1, BMI = 1, educat = 1, 
                                leisurePA = 1, placebirth = 1, diabd = 1, ahi = 1, 
                                odi = 1, t90 = 1, shannon = 1, diabd = 1, hypertension = 1, 
                                Fibrer = 1, Energi_kcal = 0, ppi = 1, hypermed = 1, 
                                dyslipmed = 1))
  
  export2xls(t1, file = paste0(descriptive.folder,'Table1.xlsx'))
  
  wrf <- "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/"
  
  export2xls(t1, file = paste0(wrf,'Table1.xlsx'))
  
  
  # T90 groups 
  t = compareGroups(t90cat ~ age + Sex + smokestatus + Alkohol + BMI + educat +
                      leisurePA + placebirth + diabd + hypertension + dyslipidemia + metformin + 
                      hypermed + dyslipmed+ ppi+Fibrer+Energi_kcal+ahi+ odi+ t90+ shannon, 
                    data = datatable1[datatable1$valid.t90=="yes",], 
                    include.miss = FALSE, chisq.test.perm = TRUE, 
                    method = c(age = 2, Alkohol = 2, BMI = 2, 
                               ahi = 2 , odi = 2, t90 = 2, 
                               shannon = 2, Fibrer = 2, Energi_kcal = 2))
  
  t1 <-  createTable(t, hide.no = "no", show.p.overall = FALSE, show.all=T, 
                     digits = c(age = 1, Sex = 1, Alkohol = 1, BMI = 1, educat = 1, 
                                leisurePA = 1, placebirth = 1, diabd = 1, ahi = 1, 
                                odi = 1, t90 = 1, shannon = 1, diabd = 1, hypertension = 1, 
                                Fibrer = 1, Energi_kcal = 0, ppi = 1, hypermed = 1, 
                                dyslipmed = 1))
  
  export2xls(t1, file = paste0(descriptive.folder,'Table1_t90cat.xlsx'))
  
  wrf <- "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/"
  
  export2xls(t1, file = paste0(wrf,'Table1_t90cat.xlsx'))
  
  
  # ODI groups 
  t = compareGroups(odicat ~ age + Sex + smokestatus + Alkohol + BMI + educat +
                      leisurePA + placebirth + diabd + hypertension + dyslipidemia + metformin + 
                      hypermed + dyslipmed+ ppi+Fibrer+Energi_kcal+ahi+ odi+ t90+ shannon, 
                    data = datatable1[datatable1$valid.t90=="yes",], 
                    include.miss = FALSE, chisq.test.perm = TRUE, 
                    method = c(age = 2, Alkohol = 2, BMI = 2, 
                               ahi = 2 , odi = 2, t90 = 2, 
                               shannon = 2, Fibrer = 2, Energi_kcal = 2))
  
  t1 <-  createTable(t, hide.no = "no", show.p.overall = FALSE, show.all=T, 
                     digits = c(age = 1, Sex = 1, Alkohol = 1, BMI = 1, educat = 1, 
                                leisurePA = 1, placebirth = 1, diabd = 1, ahi = 1, 
                                odi = 1, t90 = 1, shannon = 1, diabd = 1, hypertension = 1, 
                                Fibrer = 1, Energi_kcal = 0, ppi = 1, hypermed = 1, 
                                dyslipmed = 1))
  
  export2xls(t1, file = paste0(descriptive.folder,'Table1_odicat.xlsx'))
  
  wrf <- "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/"
  
  export2xls(t1, file = paste0(wrf,'Table1_odicat.xlsx'))

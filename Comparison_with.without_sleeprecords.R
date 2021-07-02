# Project: Sleep apnea and microbiota

# Comparing those with sleep records to those without

# Loading packages
pacman::p_load(tidyverse, Hmisc, sjmisc, compareGroups, data.table)

rm(list=ls())

# Import sleep recording data 
sleep = fread("/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording/Data_Sleep_SCAPIS_Uppsala.csv", 
              header=T, na.strings=c("", "NA"))
sleep = sleep[flutv4h==1 & o2utv4h==1,]

# Uploading phenotype data ####
pheno=fread("/home/baldanzi/Datasets/sleep_SCAPIS/SCAPIS-DATA-PETITION-170-20210315.csv",
            header = T, na.strings=c("", "NA"))
# Restricting dataset to Uppsala participants only 
pheno = pheno[Site=="Site5", ]

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
  
  # Self-reported Sleep apnea diagnosis 
  pheno$apnea_self = factor(pheno$cqhe058,  levels = c("NO", "YES"), labels = c("no", "yes"))

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

  
  # Metabolon data ####
  metabolon=fread("/home/baldanzi/Datasets/Metabolon/clean/scapis_metabolon_scapisid.tsv", header=T)
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
  
  
# Adjusting ids 
  setnames(sleep, "id", "SCAPISid") #renaming the id variable
  setnames(pheno, "Subject", "SCAPISid") #renaming the id variable
  
# Variable for identifying those with sleep recording
  pheno[, sleep.recording:="no"]
  pheno[SCAPISid %in% sleep$SCAPISid, sleep.recording:="yes"]
  
# Merge medication ####
  pheno <- merge(pheno, medication, by="SCAPISid", all.x=T)
  
   #Compare groups   
  label(pheno[["educat"]]) <- "Education"
  label(pheno[["leisurePA"]]) <- "Leisure Physical Act."
  label(pheno[["ppi"]]) <- "PPI"
  label(pheno[["smokestatus"]]) <- "Smoking status"
  label(pheno[["Alkohol"]]) <- "Alcohol (g/month)"
  label(pheno[["diabd"]]) <- "Diabetes"
  label(pheno[["hypertension"]]) <- "Hypertension"
  label(pheno[["hypermed"]]) <- "Antihypertensive med."
  label(pheno[["Fibrer"]]) <- "Fiber"
  label(pheno[["Energi_kcal"]]) <- "Total energy intake"
  label(pheno[["dyslipidemia"]]) <- "Dyslipidemia"
  label(pheno[["apnea_self"]]) <- "Apnea diagnosis"
  
  t = compareGroups(sleep.recording ~ age + Sex + smokestatus + Alkohol + BMI + educat +
                      leisurePA + diabd + hypertension + dyslipidemia + metformin + 
                      hypermed + dyslipmed+ ppi+Fibrer+Energi_kcal+apnea_self, data= pheno, 
                    include.miss = FALSE, chisq.test.perm = TRUE)
  t1 = createTable(t, hide.no = "no")
  t1
  export2html(t1, file='/home/baldanzi/Sleep_apnea/Descriptive/with.without_sleeprec.html',
              header.labels = c(p.overall = "p-value"))
  
  
  
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Script for data management and merging the sleep and MGS datasets

# Data management ####

# Loading packages
pacman::p_load(tidyverse, grid, chron, rio, Hmisc, sjmisc, summarytools, data.table)

rm(list=ls())

  # Uploading phenotype data ####
  pheno=fread("/home/baldanzi/Datasets/sleep_SCAPIS/SCAPIS-DATA-PETITION-170-20210315.csv", header = T, na.strings=c("", "NA"))
  # Restricting dataset to Uppsala participants only 
  pheno = pheno[Site=="Site5", ]


  # Uploading MGS data ####
  data.MGS = fread("/home/baldanzi/Datasets/MGS/clean/MGS_relative_abundance_4839_upp_4980_malmo.tsv")

  # Restricting to Uppsala participants 
  data.MGS = data.MGS[SCAPISid %in% grep("5-",data.MGS$SCAPISid, value = T ),] 
  # 4839 Uppsala participants have gut microbiome data
  
  # Fixing MGS names to latest annotation 
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  
  mgs.names.index <- grep("____",names(data.MGS))
  names(data.MGS)[mgs.names.index] <- taxonomy$maintax_mgs
  

  pob = fread("/home/baldanzi/Datasets/Placeofbirth.csv", header=T, sep=",")
  pob$placebirth = factor(pob$q005a, 
                        levels = c("scandinavia", "europe", "asia", "other"), 
                        labels = c("Scandinavia", "Europe", "Asia", "other"))
  pob[,q005a:=NULL]
  #restricting to Uppsala 
  pob = pob[pob$SCAPISid %in% grep("5-",pob$SCAPISid, value = T ),]

  # Metabolon data ####
  metabolon=fread("/home/baldanzi/Datasets/Metabolon/clean/scapis_metabolon_batchnorm_scapisid.tsv", header=T)

  # Metabolon Metadata - collection date ####
  # Import metabolon metadata 

  metab_metadata=fread('/home/baldanzi/Datasets/Metabolon/clean/metabolon_metadata_scapisid.tsv', header=T)
  a = c("SCAPISid", "COLLECTION_DATE")
  metab_collection_date = metab_metadata[,a,with=F]
  setnames(metab_collection_date,"COLLECTION_DATE", "metabolon_collection_date")

  # 28 participants have sleep and gut microbiota data but do not have metabolomics data


  # Recoding variables ####

  #Month of anthropometric collection date
  pheno[,visit.month:=format(as.POSIXct(pheno$AnthropometryCollectionDate),"%B")]
  pheno[,visit.month:=factor(visit.month,c(month.name,"June.July"))]
  
      # Merging June to July due to the low number of July participants 
      #pheno[visit.month %in% c("June","July"),visit.month:="June.July"]
      
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

  # remove extremes outlieers 
  pheno[Energi_kcal<500 | Energi_kcal >6000, Energi_kcal:=NA] # Removing invalid values
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


# Medication data based on metabolomics measurement
a = c("SCAPISid","MET_100002725","MET_100002808","MET_100002405")
#b = c("SCAPIS_RID", "MET_100002725","MET_100002808","MET_100002405")
medication = metabolon[,a,with=F]
medication[,ppi:=ifelse(MET_100002725>min(MET_100002725,na.rm=T) |
                      MET_100002808>min(MET_100002808,na.rm=T),1,0)]
medication[,metformin:=ifelse(MET_100002405>min(MET_100002405,na.rm=T),1, 0)]


medication[,ppi:=factor(ppi, levels = c(0,1), 
                 labels = c("no", "yes"))]

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


# Merging pheno (phenotype data) with the remaining data
nrow(pheno) # 5036 Uppsala individuals with phenotype data
setnames(pheno, "Subject", "SCAPISid") #renaming the id variable
nrow(data.MGS) # 4839 Uppsala Participants 

#  pheno + data.MGS
sum(pheno$SCAPISid %in% data.MGS$SCAPISid) # 4839 with pheno and MGS data 
pheno = merge(pheno, data.MGS, by = "SCAPISid", all.x=F, all.y=T)

#  pheno + data.MGS + place of birth
pheno = merge(pheno, pob, by = "SCAPISid", all.x=T, all.y=F)

# pheno + data.MGS + place of birth + diet 
#pheno = merge(pheno, diet, by ="SCAPISid", all.x=T, all.y=F)

# ... + medication 
pheno = merge(pheno,medication, by="SCAPISid", all.x=T, all.y=F)

# ... + metabolon collection date
pheno=merge(pheno, metab_collection_date, by="SCAPISid", all.x=T, all.y=F)

  # 51 participants have pheno and gut microbiota data but do not have metabolomics data
  sum(pheno[,SCAPISid] %in% metabolon[,SCAPISid])
  nrow(pheno)-sum(pheno[,SCAPISid] %in% metabolon[,SCAPISid])


  # Changing order
  setcolorder(pheno, c("SCAPISid", "age", "Sex"))

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

  # merging CM variables with the pheno. 
  a = c("SCAPISid", "plate", "received", "read.pairs.per.sample", "dna.yield")
  pheno = merge(pheno, cmvar[,a,with=F], by="SCAPISid", all=F, all.x=T, all.y=F)
  pheno$plate = as.factor(pheno$plate)
  pheno$received = as.factor(pheno$received)

  # Import sleep recording data 
  sleep <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording/sleep.rds")
  
  pheno <- merge(pheno, sleep, by='SCAPISid', all.x=T, all.y=F)
  
  # Exclude CPAP users 
  pheno[cqhe061=='CPAP', c("valid.ahi","valid.t90"):="no"]
  pheno[cpap=='yes', c("valid.ahi","valid.t90"):="no"]
  pheno[valid.ahi=='no',ahi:=NA]
  pheno[valid.t90=='no',t90:=NA]
  
  # T90 categories 
  pheno[t90!=0, t90cat :=  cut(t90,breaks = quantile(t90,
                                                     probs = seq(0,1,by=.33),
                                                     na.rm=T), include.lowest = T)]
  pheno[,t90cat:=factor(t90cat, levels = c("0", levels(t90cat)))]
  pheno[t90==0, t90cat := '0' ]

  # Save pheno data ####
  saveRDS(pheno, file="/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

  #valid.ahi + pheno
  nrow(sleep[valid.ahi=='yes',]) # 3301 individuals with valid flow and sat monitoring 
  sum(sleep[valid.ahi=='yes',]$SCAPISid %in% pheno$SCAPISid) # 3208 with valid monitoring and pheno and MGS data 
  valid.ahi = merge(sleep[valid.ahi=='yes',], pheno, by = "SCAPISid", all=F )

  #Saving the data #### 
  write.table(valid.ahi, file = "/home/baldanzi/Datasets/sleep_SCAPIS/validsleep.MGS.Upp.tsv",row.names = FALSE,
              sep = '\t',quote=FALSE)

  # valid.t90 + pheno
  valid.t90 = merge(sleep[valid.t90=='yes',], pheno, by = "SCAPISid", all=F )

  # Save valid.t90 data
  write.table(valid.t90, file = "/home/baldanzi/Datasets/sleep_SCAPIS/validodi.MGS.Upp.tsv", row.names = FALSE,
              sep = "\t",quote = FALSE)

  # 3622 individuals have a valid T90 measurement 




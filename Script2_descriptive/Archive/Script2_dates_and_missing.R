#---------------------------------------------------------------------------#
# We wanted to check if blood samples for metabolomics and Sleep records (ApneaLink)
#took place close to each other (difference in calender date between the two)


#### Dates ####
# Was sleep apnea assessed close to when metabolomics was collected? 
valid.ahi[,metabolon_collection_date := as.Date(metabolon_collection_date, 
                                                format = "%d-%b-%y")]

# range and mean difference between the two dates 
diff=abs(valid.ahi[,Date]-valid.ahi[,metabolon_collection_date])
summary(as.numeric(diff))
sum(diff>30,na.rm=T) # 0 individuals have a difference greater than 
# 30 days between dates

# Sleep assessment dates and Antropometric Collection Date 
# Missing data in anthropometric collection 
sum(is.na(valid.ahi[,AnthropometryCollectionDate])) # ZERO
# difference between the two dates 
diff2=abs(as.Date(valid.ahi[,Date])-as.Date(valid.ahi[,AnthropometryCollectionDate]))
summary(as.numeric(diff2))
sum(diff2>30,na.rm=T) # 23 had a difference greater than 30 days

#--------------------------------------------------------------------------#
#### Missingness #### 
message("")
message("Table of complete vs incomplete observations")
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

saveRDS(t1, file='/home/baldanzi/Sleep_apnea/Descriptive/completVSincompletObs.rds')


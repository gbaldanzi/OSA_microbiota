# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# Last update: 2021-10-04

  #Variables are: SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, WaistHip, educat,
#leisurePA, placebirth, diabd, hypertension, dyslipidemia, diabmed, 
#hypermed, dyslipmed, ppi, fiber, EI,ESS,apnea_self, apneatto_self,
#cpap_self, splint_self, apneasurgery_self,
#ahi, odi, sat90, cpap, splint 
  
dat1 = pheno[valid.ahi=='yes',] %>% select(SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, educat,
                                           leisurePA, placebirth, diabd, hypertension, dyslipidemia, metformin, 
                                           hypermed, dyslipmed, ppi, Fibrer, Energi_kcal, ahi, t90, shannon)


#dat1 = pheno %>% select(SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, educat,
#                                           leisurePA, placebirth, diabd, hypertension, dyslipidemia, metformin, 
 #                                          hypermed, dyslipmed, ppi, Fibrer, Energi_kcal, ahi, t90, shannon)


  dat1[OSAcat == "no OSA", OSAcat:= "No Sleep Apnea"]
  dat1[,OSAcat:=factor(OSAcat, levels = c("No Sleep Apnea", "Mild", "Moderate", "Severe"))]

# Variable diagnostics ####
# Number of individuals 
nrow(dat1) #3175



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
             "Anti-hypertensive med.", "Cholesterol-lowering med.", "PPI","Fiber (g/day)","Energy intake (kcal/day)",
             "AHI (events/h)", "T90 (%)", "Shannon Index")

j=1
for(i in names(datatable1)){
  label(datatable1[[i]]) <- mylabel[j]
  j=j+1
}

# Create table of population characteristics by OSA group (no OSA, Mild, Moderate, or Severe OSA)
t = compareGroups(OSAcat ~ age + Sex + smokestatus + Alkohol + BMI + educat +
                    leisurePA + placebirth + diabd + hypertension + dyslipidemia + metformin + 
                    hypermed + dyslipmed+ ppi+Fibrer+Energi_kcal+ahi+ t90+ shannon, 
                  data= datatable1, 
                  include.miss = FALSE, chisq.test.perm = TRUE)
t1 = createTable(t, hide.no = "no")

saveRDS(t1, file='/home/baldanzi/Sleep_apnea/Descriptive/sleepapnea_table1.rds')

export2xls(t1, file='/home/baldanzi/Sleep_apnea/Descriptive/Table1.xlsx')

export2xls(t1, file='Table1.xlsx')


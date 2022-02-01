# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# Last update: 2022-02-01

# Script to create a Table 1 with participants baseline characteristics 

  # Folder where the table will be outputed 
  output.plot="/home/baldanzi/Sleep_apnea/Descriptive/"

# Table 1 only includes those participants with valid AHI 

    #Variables are: SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, WaistHip, educat,
  #leisurePA, placebirth, diabd, hypertension, dyslipidemia, diabmed, 
  #hypermed, dyslipmed, ppi, fiber, EI,ESS,apnea_self, apneatto_self,
  #cpap_self, splint_self, apneasurgery_self,
  #ahi, odi, sat90, cpap, splint 
  
  dat1 = pheno[valid.ahi=='yes',] %>% select(SCAPISid, OSAcat , age, Sex, smokestatus, Alkohol, BMI, educat,
                                           leisurePA, placebirth, diabd, hypertension, dyslipidemia, metformin, 
                                           hypermed, dyslipmed, ppi, Fibrer, Energi_kcal, ahi, odi,t90, shannon)


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
             "AHI (events/h)","ODI (events/h)", "T90 (%)", "Shannon Index")

  j=1
  for(i in names(datatable1)){
    label(datatable1[[i]]) <- mylabel[j]
    j=j+1
  }

  # Create table of population characteristics by OSA group (no OSA, Mild, Moderate, or Severe OSA)
  t = compareGroups(OSAcat ~ age + Sex + smokestatus + Alkohol + BMI + educat +
                    leisurePA + placebirth + diabd + hypertension + dyslipidemia + metformin + 
                    hypermed + dyslipmed+ ppi+Fibrer+Energi_kcal+ahi+ odi+ t90+ shannon, 
                  data= datatable1, 
                  include.miss = FALSE, chisq.test.perm = TRUE)
  t1 = createTable(t, hide.no = "no")
  

  saveRDS(t1, file = paste0(output.plot,'sleepapnea_table1.rds'))
  

  export2xls(t1, file = paste0(output.plot,'Table1.xlsx'))

# Script 5 - Main 

# Script to run the permutation analysis between OSA/AHI/T90 and BC/AD

# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)


#Covariates 
  # basic.model : adjust for age + sex + alcohol + smoking + plate + received 
  basic.model<-   c("age", "Sex", "Alkohol","smokestatus","plate")

  # full.model = basic.model + fiber intake + Energy intake + physical activity + education + country of birth + ppi + metformin +  antihypertensive + cholesterol-lowering 
  full.model <-  c(basic.model, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth", "visit.month", "metformin","hypermed","dyslipmed","ppi")
  
  # full.model_BMI = full.model + BMI
  full.model_BMI <-  c(full.model, "BMI")

  
  #source('Script5_permanova/perma_ahi_bc/perma_ahi_bc_main.R')

  source('Script5_permanova/perma_odi_bc/perma_odi_bc_main.R')

  #source('Script5_permanova/perma_t90_bc/perma_t90_bc_main.R')


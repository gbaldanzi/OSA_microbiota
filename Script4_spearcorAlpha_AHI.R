# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# This code will investigate the alpha diversity in relation AHI, OSAcat, and T90%.  
# Three models will be used 
# Results will be saved in three files depending on the exposure (AHI,OSAcat,T90%)


# Loading packages 
pacman::p_load(data.table,ppcor, fastDummies, vegan, parallel)

# Cleaning the environment 
rm(list = ls())

# Defining output folders 
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.ahi = fread("validsleep_MGS.shannon.BC_Upp.tsv",header=T, na.strings=c("", "NA")) 
setnames(valid.ahi, "pob", "placebirth")

#Spearman correlation function ####
spearman.function = function(x1, x2, covari=NULL, data){
  if(any(class(data) %in% "data.table")){setDF(data)}
  numeric.covari = names(data[covari])[sapply(data[covari],is.numeric)]
  factor.covari =  names(data[covari])[sapply(data[covari],is.factor)]
  factor.covari = c(factor.covari,
                    names(data[covari])[sapply(data[covari],is.character)])
  #Exclude rows with incomplete observation
  if(any(is.na(data[,x1]))){stop("NA in x1")}
  if(any(is.na(data[,x2]))){stop("NA in x2")}
  
  notexclude = which(apply(data[c(numeric.covari,factor.covari)], 1, function(x){all(!is.na(x))}))
  temp.data=data[notexclude,]
  
  if(length(factor.covari)>0){  #Dummies variables 
    temp.data = dummy_cols(temp.data, select_columns = factor.covari,
                           remove_most_frequent_dummy = T, 
                           remove_selected_columns = T)
    factor.covari = names(temp.data)[!names(temp.data) %in% names(data)]
    #covariates for final model 
    cov=c(numeric.covari,factor.covari)
  }
  #Partial Sperman correlation 
  result=data.frame(matrix(ncol=7, nrow=length(x1)))
  for(i in 1:length(x1)){
    res=pcor.test(temp.data[,x1[i]],temp.data[,x2],temp.data[,cov],
                  method = "spearman")
    a=data.frame(x=x1[i],y=x2,
                 estimate=res$estimate,
                 p.value=res$p.value,
                 N=res$n,
                 method = res$Method,
                 covariates=paste(covari,collapse = "+"))
    result[i,]=a
  }
  if(length(x1)>1){
    result$q.value = p.adjust(result[,4], method="BH")
  }
  return(result)
}

#-----------------------------------------------------------------------------#

# Correlation between AHI and Shannon ####

# Transforming two-level factor variables into numeric variables 
dades = copy(valid.ahi)
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming factor variables 
dades[,plate:=as.factor(dades$plate)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

# Preparing exposure, outcomes, and covariates
exposure="ahi"
outcomes="shannon"

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate")
# model 2 = model 1 + BMI 
model2 <-  c(model1,"BMI")
# model 3 = model 2 + fiber intake +energy intake+ physical activity + education + country of birth 
model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth")
# OLD model 3 = model 2 + diabetes + hypertension + dyslipidemia, medication 
# OLD model3 = c(model2,"diabd","hypertension","dyslipidemia","diabmed","hypermed","dyslipmed","ppi")

listmodels=list(model1,model2,model3)

# Run Spearman correlation for the models.
res.alpha = lapply(listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)})
names(res.alpha) = c("model1", "model2", "model3")

res.model1 = res.alpha[[1]]
res.model2 = res.alpha[[2]]
res.model3 = res.alpha[[3]]

res.model1$model= "model1"
res.model2$model= "model2"
res.model3$model= "model3"

res.alpha = as.data.table(rbind(res.model1, res.model2, res.model3))

names(res.alpha) = c("MGS", "exposure", "cor.coeficient", "p.value", 
                     "N", "method", "covariates","model")
#----------------------------------------------------------------------------#
# Sensitivity analysis - remove individuals who use medication 
dades = copy(valid.ahi)
dades <-  dades[dades$ppi == "no",] #60
dades <-  dades[dades$metformin == "no",] #57 
dades <-  dades[dades$hypermed == "no",] #593
dades <-  dades[dades$dyslipmed == "no",] # 239
nrow(dades) #2318

# Transforming factor variables 
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
dades[,plate:=as.factor(dades$plate)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

# Spearman correlation 
res.sa = spearman.function(x1=outcomes,x2=exposure,covari = model3,data = dades)

res.sa$model= "model3_noMedication"

#naming coluns
names(res.sa) <- c("MGS", "exposure", "cor.coeficient", "p.value", 
                "N", "method", "covariates","model")

res.alpha = rbind(res.alpha,res.sa)

fwrite(res.alpha, file = paste0(output,"cor_ahi_alpha.tsv"), sep="\t")

#-----------------------------------------------------------------------------#

# Correlation between AHI and Shannon by BMI group ####

# the correlation between AHI and shannon will be separetely by the 3 BMI groups:
# <25, >=25 & <30, and >=30.

# Transforming two-level factor variables into numeric variables 
dades = copy(valid.ahi)
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming factor variables 
dades[,plate:=as.factor(dades$plate)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

# Preparing exposure, outcomes, and covariates
exposure="ahi"
outcomes="shannon"

# By BMI group 
for(group in unique(valid.ahi[,BMIcat])){
  dades2 = dades[BMIcat==group,]


# Run Spearman correlation for the models.
res.alpha = lapply(listmodels,function(mod){
     spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades2)})
names(res.alpha) = c("model1", "model2", "model3")
  
res.model1 = res.alpha[[1]]
res.model2 = res.alpha[[2]]
res.model3 = res.alpha[[3]]



res.model1$model= "model1"
res.model2$model= "model2"
res.model3$model= "model3"

res.model1$bmi= group
res.model2$bmi= group
res.model3$bmi= group

res.alpha = as.data.table(rbind(res.model1, res.model2, res.model3))

names(res.alpha) = c("MGS", "exposure", "cor.coeficient", "p.value", 
                     "N", "method", "covariates","model","bmi")
print(res.alpha)

#----------------------------------------------------------------------------#
print("Sensitivity analysis")
# Sensitivity analysis - remove individuals who use medication 
dades = copy(valid.ahi)
dades <-  dades[dades$ppi == "no",] #60
dades <-  dades[dades$metformin == "no",] #57 
dades <-  dades[dades$hypermed == "no",] #593
dades <-  dades[dades$dyslipmed == "no",] # 239
nrow(dades) #2318

# Transforming factor variables 
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
dades[,plate:=as.factor(dades$plate)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

# BMI group
dades2 <- dades[BMIcat==group,]

# Spearman correlation 
res.sa = spearman.function(x1=outcomes,x2=exposure,covari = model3,data = dades2)

res.sa$model= "model3_noMedication"
res.sa$bmi= group

#naming coluns
names(res.sa) <- c("MGS", "exposure", "cor.coeficient", "p.value", 
                   "N", "method", "covariates","model","bmi")

res.alpha = rbind(res.alpha,res.sa)

fwrite(res.alpha, file = paste0(output,"cor_ahi_alpha_bmi",group,".tsv"), sep="\t")
}
#-----------------------------------------------------------------------------#

# Correlation between T90% and Shannon ####

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.odi = fread('/home/baldanzi/Datasets/sleep_SCAPIS/validodi.MGS.Upp.tsv',
                  header=T, na.strings=c("", "NA")) 
setnames(valid.odi, "pob", "placebirth")
setnames(valid.odi, "sat90", "T90")
valid.odi = valid.odi[!is.na(T90),]


# Calculating alpha-diversity 
a = grep("____",names(valid.odi),value=T) # vector with MGS names 
valid.odi$shannon=diversity(valid.odi[,a, with=F],index="shannon") #estimating shannon per individual



# Transforming two-level factor variables into numeric variables 
dades = copy(valid.odi)
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming factor variables 
dades[,plate:=as.factor(dades$plate)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

# Preparing exposure, outcomes, and covariates
exposure="T90"
outcomes="shannon"

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate")
# model 2 = model 1 + BMI 
model2 <-  c(model1,"BMI")
# model 3 = model 2 + fiber intake +energi intake+ physical activity + education + country of birth 
model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth")
# OLD model 3 = model 2 + diabetes + hypertension + dyslipidemia, medication 
# OLD model3 = c(model2,"diabd","hypertension","dyslipidemia","diabmed","hypermed","dyslipmed","ppi")

listmodels=list(model1,model2,model3)

# Run Spearman correlation for the models.
res.alpha = lapply(listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)})
names(res.alpha) = c("model1", "model2", "model3")

res.model1 = res.alpha[[1]]
res.model2 = res.alpha[[2]]
res.model3 = res.alpha[[3]]

res.model1$model= "model1"
res.model2$model= "model2"
res.model3$model= "model3"

res.alpha = as.data.table(rbind(res.model1, res.model2, res.model3))

names(res.alpha) = c("MGS", "exposure", "cor.coeficient", "p.value", 
                     "N", "method", "covariates","model")
#----------------------------------------------------------------------------#
# Sensitivity analysis - remove individuals who use medication 
dades = copy(valid.odi)
dades <-  dades[dades$ppi == "no",] #60
dades <-  dades[dades$metformin == "no",] #57 
dades <-  dades[dades$hypermed == "no",] #593
dades <-  dades[dades$dyslipmed == "no",] # 239
nrow(dades) #2318

# Transforming factor variables 
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
dades[,plate:=as.factor(dades$plate)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

# Spearman correlation 
res.sa = spearman.function(x1=outcomes,x2=exposure,covari = model3,data = dades)

res.sa$model= "model3_noMedication"

#naming coluns
names(res.sa) <- c("MGS", "exposure", "cor.coeficient", "p.value", 
                   "N", "method", "covariates","model")

res.alpha = rbind(res.alpha,res.sa)

#Saving results 
fwrite(res.alpha, file = paste0(output,"cor_t90_alpha.tsv"), sep="\t")

#-----------------------------------------------------------------------------#

# Correlation between T90% and Shannon by BMI group ####

# the correlation between AHI and shannon will be separetely by the 3 BMI groups:
# <25, >=25 & <30, and >=30.

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.odi = fread('/home/baldanzi/Datasets/sleep_SCAPIS/validodi.MGS.Upp.tsv',
                  header=T, na.strings=c("", "NA")) 
setnames(valid.odi, "pob", "placebirth")
setnames(valid.odi, "sat90", "T90")
valid.odi = valid.odi[!is.na(T90),]


# Calculating alpha-diversity 
a = grep("____",names(valid.odi),value=T) # vector with MGS names 
valid.odi$shannon=diversity(valid.odi[,a, with=F],index="shannon") #estimating shannon per individual

# Transforming two-level factor variables into numeric variables 
dades = copy(valid.odi)
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming factor variables 
dades[,plate:=as.factor(dades$plate)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

# Preparing exposure, outcomes, and covariates
exposure="T90"
outcomes="shannon"

#By BMI group
for(group in unique(valid.ahi[,BMIcat])){

dades2 = dades[BMIcat==group,]
  
# Run Spearman correlation for the models.
res.alpha = lapply(listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades2)})
names(res.alpha) = c("model1", "model2", "model3")

res.model1 = res.alpha[[1]]
res.model2 = res.alpha[[2]]
res.model3 = res.alpha[[3]]

res.model1$model= "model1"
res.model2$model= "model2"
res.model3$model= "model3"

res.model1$bmi= group
res.model2$bmi= group
res.model3$bmi= group

res.alpha = as.data.table(rbind(res.model1, res.model2, res.model3))

names(res.alpha) = c("MGS", "exposure", "cor.coeficient", "p.value", 
                     "N", "method", "covariates","model","bmi")
#----------------------------------------------------------------------------#
print("Sensitivity analysis (no medication users) correlation T90% and Shannon Index")
# Sensitivity analysis - remove individuals who use medication 
dades = copy(valid.odi)
dades <-  dades[dades$ppi == "no",] #60
dades <-  dades[dades$metformin == "no",] #57 
dades <-  dades[dades$hypermed == "no",] #593
dades <-  dades[dades$dyslipmed == "no",] # 239
nrow(dades) #2318

# Transforming factor variables 
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
dades[,plate:=as.factor(dades$plate)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

dades2 = dades[BMIcat==group,]

# Spearman correlation 
res.sa = spearman.function(x1=outcomes,x2=exposure,covari = model3,data = dades2)

res.sa$model= "model3_noMedication"
res.sa$bmi= group

#naming coluns
names(res.sa) <- c("MGS", "exposure", "cor.coeficient", "p.value", 
                   "N", "method", "covariates","model","bmi")

res.alpha = rbind(res.alpha,res.sa)

#Saving results 
fwrite(res.alpha, file = paste0(output,"cor_t90_alpha_bmi",group,".tsv"), sep="\t")
}
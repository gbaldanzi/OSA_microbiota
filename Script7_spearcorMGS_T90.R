# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# This code will investigate the association between MGS and T90% using 3 different models. 
# Results are save in one single .tsv file with a column "model" identifying the models. 
# A important part of this code is a Partial Spearman correlation to be run in 
#parallel nodes (1 for each model). Therefore, this code requires at least 3 nodes. 

# Results saved at "/home/baldanzi/Sleep_apnea/Results/cor_t90_mgs.tsv"

rm(list = ls())
# Loading packages 
pacman::p_load(data.table,ppcor, fastDummies, vegan,parallel)

# Output folders 
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.odi = fread("validodi.MGS.Upp.tsv",header=T, na.strings=c("", "NA")) 
setnames(valid.odi, "pob", "placebirth")
setnames(valid.odi, "sat90", "T90")
valid.odi = valid.odi[!is.na(T90),]

# Calculating alpha-diversity 
a = grep("____",names(valid.odi),value=T) # vector with MGS names 
valid.odi$shannon=diversity(valid.odi[,a, with=F],index="shannon") #estimating shannon per individual

# Transforming two-level factor variables into numeric variables 
dades = copy(valid.odi)
#a= c("Sex","ppi", "diabd","hypertension","dyslipidemia","diabmed","hypermed","dyslipmed")
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming factor variables 
dades[,plate:=as.factor(dades$plate)]
dades[,received:=as.factor(dades$received)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

#Spearman correlation function ####
spearman.function = function(x1, x2, covari=NULL, data){
  if(any(class(data) %in% "data.table")){setDF(data)}
  numeric.covari = names(data[covari])[sapply(data[covari],is.numeric)]
  factor.covari =  names(data[covari])[sapply(data[covari],is.factor)]
  factor.covari = c(factor.covari,
                    names(data[covari])[sapply(data[covari],is.character)])
  #Exclude rows with incomplete observation
  if(any(c(is.na(data[,x1]),is.na(data[,x2])))){print("NA in x1 or x2");stop()}
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

# Correlation between T90 and MGS - Step1 ####

#Preparing exposure, outcomes, and covariates
exposure="T90"
outcomes=grep("___",names(valid.odi),value=T)

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate + received 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","received","shannon")
# model 2 = model 1 + BMI 
model2 <-  c(model1,"BMI")
# model 3 = model 2 + fiber intake + physical activity + education + country of birth 
model3 <-  c(model2, "Fibrer", "leisurePA", "educat","placebirth")
# OLD model 3 = model 2 + diabetes + hypertension + dyslipidemia, medication 
# OLD model3 = c(model2,"diabd","hypertension","dyslipidemia","diabmed","hypermed","dyslipmed","ppi")

listmodels=list(model1,model2,model3)

# Preparing parallelism
c1 = makeCluster(3)
clusterEvalQ(c1, library(vegan))
clusterEvalQ(c1, library(data.table))
clusterEvalQ(c1, library(ppcor))
clusterEvalQ(c1, library(fastDummies))
clusterExport(c1, c("exposure","outcomes", "dades", "model1", "model2", "model3","spearman.function"))
step1 = parLapply(c1, listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)})
names(step1) = c("model1", "model2", "model3")
step1 = lapply(step1,function(x){x[order(x$q.value),]})

res.model1 = step1[[1]]
res.model2 = step1[[2]]
res.model3 = step1[[3]]

res.model1$model= "model1"
res.model2$model= "model2"
res.model3$model= "model3"

step1.res = as.data.table(rbind(res.model1, res.model2, res.model3))
stopCluster(c1)

names(step1.res) = c("MGS", "exposure", "cor.coeficient", "p.value", 
                     "N", "method", "covariates","q.value","model")


# Sensitivity analysis - remove individuals who use medication 
dades <-  copy(valid.odi)
dades <-  dades[ppi == "no",] #60
dades <-  dades[metformin == "no",] #57 
dades <-  dades[hypermed == "no",] #593
dades <-  dades[dyslipmed == "no",] # 239
nrow(dades) #2318

# Transforming two-level factor variables into numeric variables 
a= c("Sex")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming factor variables 
dades[,plate:=as.factor(dades$plate)]
dades[,received:=as.factor(dades$received)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

#Prepare data.frame data will receive the results 
res = data.frame(matrix(ncol=7, nrow=1985))

# Spearman correlation 
res = spearman.function(x1=outcomes,x2=exposure,covari = model3,data = dades)

res$model= "model3_noMedication"

# Sort by q.value 
res <- res[order(res$q.value),]
#naming coluns
names(res) <- c("MGS", "exposure", "cor.coeficient", "p.value", 
                "N", "method", "covariates","q.value","model")

step1.res = rbind(step1.res,res)

fwrite(step1.res, file = paste0(output,"cor_t90_mgs.tsv"), sep="\t")


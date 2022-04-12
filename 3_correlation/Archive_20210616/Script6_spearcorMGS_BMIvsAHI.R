# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# This code aims to find discordant MGS associated with either BMI or AHI


# Results saved at "/home/baldanzi/Sleep_apnea/Results/discordant_mgs.tsv"

  rm(list = ls())
# Loading packages 
  pacman::p_load(data.table,ppcor, fastDummies, vegan, ggplot2,parallel)

# Output folders 
  output = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Import data
  setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
  valid.ahi <- readRDS("validsleep_MGS.shannon_Upp.rds")
  setnames(valid.ahi, "pob", "placebirth")


#Calculating  MGS prevalence ####
  noms=grep("____",names(valid.ahi),value=T)
# presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
  data_pa <- decostand(x = valid.ahi[,noms,with=F], "pa")
# calculate sum per species
  data_sum <- data.frame(prevalence=apply(data_pa, 2, sum))
  data_sum$MGS = rownames(data_sum)
  a = data_sum$MGS[data_sum$prevalence<1] # Bacteroidales_sp.____HG3A.1976 is not present in any participant
  a = data_sum$MGS[data_sum$prevalence<5] #19 MGS are only present in less than 5 individuals
 

# Removing Bacteroidales_sp.____HG3A.1976 from the dataset 
  valid.ahi <- valid.ahi[ , -a, with=F]  # 1984 MGS remaining 

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

# Correlation between AHI and MGS NOT ADJUSTING FOR BMI ####

  # Transforming two-level factor variables into numeric variables 
  dades = copy(valid.ahi)
  a= c("Sex", "diabmed","hypermed","dyslipmed","ppi")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
  
  # Transforming factor variables 
  dades[,plate:=as.factor(dades$plate)]  
  
  # Preparing exposure, outcomes, and covariates
  exposure="ahi"
  outcomes=grep("___",names(valid.ahi),value=T)

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate + received 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")

# model 2 = model 1 + fiber intake+ Energi intake + physical activity + education + country of birth 
model2 <-  c(model1, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth")

# model 3 = model 2 + ppi + metformin + anti-hypertensive + cholesterol-lowering  
  model3 <-  c(model2, "diabmed","hypermed","dyslipmed","ppi")


listmodels=list(model1,model2, model3)

# Preparing parallelism
c1 = makeCluster(3)
clusterEvalQ(c1, library(vegan))
clusterEvalQ(c1, library(data.table))
clusterEvalQ(c1, library(ppcor))
clusterEvalQ(c1, library(fastDummies))
clusterExport(c1, c("exposure","outcomes", "dades", "model1", "model2", "model3",
                    "spearman.function"))
t0 = Sys.time()
step1 = parLapply(c1, listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)})
t1 = Sys.time()
print("Spearman correlations for AHI models 1-3 - time")
print(t1-t0)

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

  names(step1.res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")


# Sensitivity analysis - remove individuals who use medication 
  dades <-  copy(valid.ahi)
  dades <-  dades[ppi == "no",] #60
  dades <-  dades[metformin == "no",] #57 
  dades <-  dades[hypermed == "no",] #593
  dades <-  dades[dyslipmed == "no",] # 239
  nrow(dades) #2318

# Transforming two-level factor variables into numeric variables 
  a= c("Sex", "diabmed","hypermed","dyslipmed","ppi")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming factor variables 
  dades[,plate:=as.factor(dades$plate)]

#Prepare data.frame data will receive the results 
  res = data.frame(matrix(ncol=7, nrow=1985))

# Spearman correlation 
  res = spearman.function(x1=outcomes,x2=exposure,covari = model2,data = dades)

  res$model= "SA"

# Sort by q.value 
res <- res[order(res$q.value),]
#naming coluns
names(res) <- c("MGS", "exposure", "cor.coefficient", "p.value", 
                "N", "method", "covariates","q.value","model")

step1.res = rbind(step1.res,res)

fwrite(step1.res, file = paste0(output,"discordant_ahi.tsv"), sep="\t")

#--------------------------------------------------------------------------#

# Correlation between AHI and BMI with same models ####

  # Preparing Data - transform binary variables 
  dades <-  copy(valid.ahi)
  a= c("Sex", "diabmed","hypermed","dyslipmed","ppi")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

  # Transforming factor variables 
  dades[,plate:=as.factor(dades$plate)]

#Preparing exposure, outcomes, and covariates
  exposure="BMI"
  outcomes=grep("___",names(valid.ahi),value=T)

# Preparing parallelism
c1 = makeCluster(3)
clusterEvalQ(c1, library(vegan))
clusterEvalQ(c1, library(data.table))
clusterEvalQ(c1, library(ppcor))
clusterEvalQ(c1, library(fastDummies))
clusterExport(c1, c("exposure","outcomes", "dades", "model1", "model2","model3", "spearman.function"))

t0 = Sys.time()
step1 = parLapply(c1, listmodels,function(mod){
  spearman.function(x1=outcomes,x2=exposure,covari = mod,data = dades)})
t1 = Sys.time()
print("Spearman correlation BMI and MGS models 1 -2 - time")
print(t1-t0)

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

names(step1.res) = c("MGS", "exposure", "cor.coefficient", "p.value", 
                     "N", "method", "covariates","q.value","model")


# Sensitivity analysis - remove individuals who use medication 
dades <-  copy(valid.ahi)
dades <-  dades[ppi == "no",] #60
dades <-  dades[metformin == "no",] #57 
dades <-  dades[hypermed == "no",] #593
dades <-  dades[dyslipmed == "no",] # 239
nrow(dades) #2318

# Transforming two-level factor variables into numeric variables 
a= c("Sex", "diabmed","hypermed","dyslipmed","ppi")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming factor variables 
dades[,plate:=as.factor(dades$plate)]

#Prepare data.frame data will receive the results 
res = data.frame(matrix(ncol=7, nrow=1985))

# Spearman correlation 
res = spearman.function(x1=outcomes,x2=exposure,covari = model2,data = dades)

res$model= "SA"

# Sort by q.value 
res <- res[order(res$q.value),]
#naming coluns
names(res) <- c("MGS", "exposure", "cor.coefficient", "p.value", 
                "N", "method", "covariates","q.value","model")

step1.res = rbind(step1.res,res)

fwrite(step1.res, file = paste0(output,"discordant_bmi.tsv"), sep="\t")


#-----------------------------------------------------------------------------

# Merging results with taxonomy information #### 

taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
setnames(taxonomy,"maintax_mgs","MGS")

dades <- fread(paste0(output,"discordant_ahi.tsv"))
dades <- merge(dades, taxonomy, by="MGS", all.x=T)
fwrite(dades, file=paste0(output,"discordant_ahi.tsv"))

dades <- fread(paste0(output,"discordant_bmi.tsv"))
dades <- merge(dades, taxonomy, by="MGS", all.x=T)
fwrite(dades, file=paste0(output,"discordant_bmi.tsv"))
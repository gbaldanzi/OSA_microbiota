# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# This code will investigate the beta-diversity (Bray Curtis Dissimilarity and 
# Aitchison distance) in relation to OSA severity categories in 3 different models. 
# Analysis are run using a PERMANOVA approach

#parallel nodes (1 for each model). Therefore, this code requires at least 3 nodes. 

# Results saved at the folder: "/home/baldanzi/Sleep_apnea/Results/"
# File model 1 - permanova_model1_osa_ad.tsv
# File model 2 - permanova_model2_osa_ad.tsv
# File model 3 - permanova_model3_osa_ad.tsv
# File model 3 after removing medication users - permanova_model3_nomed_osa_ad.tsv

# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)

rm(list = ls())
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.ahi = fread("validsleep_MGS.shannon.BC_Upp.tsv",header=T, na.strings=c("", "NA")) 
setnames(valid.ahi, "pob", "placebirth")

# Importing AD matrix 
print("Import AD matrix")
AD <-  fread('OSA.aitchison_distmatrix.csv', header=T, sep = ',')
a <-  AD[,SCAPISid]
AD[,c("SCAPISid","OSAcat"):=NULL]
AD <- as.matrix(AD)
rownames(AD) <- a

# The PERMANOVA function 
PermanovaFunction = function(outcome, exposure, covari, data, distance = "euclidean", nodes = 1){
  require(vegan)
  require(parallel)
  require(data.table)
  if(!any(class(data) %in% "data.table")){setDT(data)}
  if(nrow(data)==0){stop("Data has nrow == 0")}
  if(ncol(data)==0){stop("Data has ncol == 0")}
  if(!any(class(get(outcome)) %in% c("matrix"))){stop("outcome should be a matrix")}
  if(length(exposure) != 1){stop("exposure should have length == 1 ")}
  
  # Permanova from vegan package cannot handle missing information 
  a = complete.cases(data[,covari,with=F])
  
  if(length(covari) > 0){
    covariates = paste(covari,collapse = "+")
    formula = as.formula(paste0(outcome,'[a,] ~',exposure,'+',covariates))
  }
  
  if(length(covari) == 0){
    formula = as.formula(paste0(outcome,'~',exposure,))
  }
  
  result <- adonis2(formula,data[a,], by="margin", method = distance ,permutations = 999, parallel = nodes)
  variables = rownames(result)
  result=as.data.frame(result)
  result[1,"N"]=sum(a)
  result$variables = variables 
  result = result[,c("variables",names(result[,-ncol(result)]))]
  return(result)
}

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

# Making sure that AD and dades have the same order of observations 
print("Match order of observations")
print(nrow(dades))
dades = dades[match(rownames(AD),dades$SCAPISid),]
print(nrow(dades))

# Outcome - character name (length=1) with matrix distance 

outc = "AD"

# Main Exposure - character name (length=1)
expo = "OSAcat"

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate + received 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate")
# model 2 = model 1 + BMI 
model2 <-  c(model1,"BMI")
# model 3 = model 2 + fiber intake + energy intake+physical activity + education + country of birth 
model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth")

# Runing PERMANOVA in parallel ####
print("prepare paralallel")
set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outc","expo","dades","model1","model2","model3"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
print("run models:1-3")
print("model1")
res1 = PermanovaFunction(outcome = outc, exposure = expo, covari = model1, data = dades, nodes = nod)
print("model2")
res2 = PermanovaFunction(outcome = outc, exposure = expo, covari = model2, data = dades, nodes = nod)
print("model3")
res3 = PermanovaFunction(outcome = outc, exposure = expo, covari = model3, data = dades, nodes = nod)

stopCluster(cl)

# Saving result
fwrite(res1, file = paste0(output,"permanova_model1_osa_ad.tsv"), sep="\t")
fwrite(res2, file = paste0(output,"permanova_model2_osa_ad.tsv"), sep="\t")
fwrite(res3, file = paste0(output,"permanova_model3_osa_ad.tsv"), sep="\t")



#---------------------------------------------------------------------------#
# Sensitivity analysis excluding medication users 
print("Sensitivity analysis")

# Sensitivity analysis - remove individuals who use medication 
dades <-  copy(valid.ahi)
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
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

# Making sure that AD and dades have the same observations 
AD = as.data.frame(AD) # Transform AD from matrix to data.frame
cols = dades[,SCAPISid] # Pass the dades observations ids to a vector
AD = AD[,cols] # Only keep the AD columns that correspond to dades observations
AD$SCAPISid = rownames(AD) # Creates a AD variable with the rownames(AD) 
AD = AD[AD$SCAPISid %in% dades[,SCAPISid], ] # Exclude AD rows that are not present in dades 
AD$SCAPISid = NULL # Exclude the SCAPISid column 

AD = as.matrix(AD) # Transform AD back to a matrix

dades = dades[match(rownames(AD),dades$SCAPISid),]  # Makes that AD and dades are in the same order

# Runing PERMANOVA in parallel ####
print("Running PERMANOVA for model with no medication")
set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outc","expo","dades","model1","model2","model3"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
res3.nomed = PermanovaFunction(outcome = outc, exposure = expo, covari = model3, data = dades, nodes = nod)

stopCluster(cl)

# Saving results 
fwrite(res3.nomed, file = paste0(output,"permanova_model3_nomed_osa_ad.tsv"), sep="\t")

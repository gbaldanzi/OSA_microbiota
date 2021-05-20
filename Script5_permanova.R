# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# This code will investigate the beta-diversity (Bray Curtis Dissimilarty and 
# Aitchison distance) in relation to AHI and T90% in 3 differents models. 
# Analysis are run using a PERMANOVA approach

#parallel nodes (1 for each model). Therefore, this code requires at least 3 nodes. 

# Results saved at "XXXXX"

# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)

rm(list = ls())
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.ahi = fread("validsleep_MGS.shannon.BC_Upp.tsv",header=T, na.strings=c("", "NA")) 
setnames(valid.ahi, "pob", "placebirth")

# Importing BC matrix 
BC = fread('OSA.BCmatrix.csv', header=T, sep = ',')
BC = as.matrix(BC)
row.names(BC) = colnames(BC)

# The PERMANOVA function 
PermanovaFunction = function(outcome, exposure, covari, data, distance = "bray", nodes = 1){
  require(vegan)
  require(parallel)
  require(data.table)
  if(!any(class(data) %in% "data.table")){setDT(data)}
  if(nrow(data)==0){stop("Data has nrow == 0")}
  if(ncol(data)==0){stop("Data has ncol == 0")}
  if(!any(class(get(outcome)) %in% c("matrix"))){stop("outcome should be a matrix ")}
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
  result[1,"N"]=length(a)
  result$variables = variables 
  result = result[,c("variables",names(result[,-ncol(result)]))]
  return(result)
}

# Transforming two-level factor variables into numeric variables 
dades = copy(valid.ahi)
a= c("Sex","ppi", "diabd","hypertension","dyslipidemia","diabmed","hypermed","dyslipmed")
dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

# Transforming factor variables 
dades[,plate:=as.factor(dades$plate)]
dades[,received:=as.factor(dades$received)]
dades[,smokestatus:=as.factor(smokestatus)]
dades[,leisurePA:=as.factor(leisurePA)]
dades[,educat:=as.factor(educat)]
dades[,placebirth:=as.factor(placebirth)]

# Making sure that BC and dades have the same order of observations 
dades = dades[match(rownames(BC),dades$SCAPISid),]

# Outcome - character name (length=1) with matrix distance 
outc = "BC"

# Main Exposure - character name (length=1)
expo = "ahi"

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate + received 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate","received")
# model 2 = model 1 + BMI 
model2 <-  c(model1,"BMI")
# model 3 = model 2 + fiber intake + physical activity + education + country of birth 
model3 <-  c(model2, "Fibrer", "leisurePA", "educat","placebirth")

# Runing PERMANOVA in parallel ####
set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outc","expo","dades","model1"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
res = PermanovaFunction(outcome = outc, exposure = expo, covari = model1, data = dades, nodes = nod)

# Saving result
fwrite(res, file = paste0(output,"permanova_model1.tsv"), sep="\t")
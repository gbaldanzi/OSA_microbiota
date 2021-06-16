# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# This code will investigate the beta-diversity (Aitchison distance) 
# in relation to T90% in 3 differents models + 1 sensitivity analysis
# Analysis are run using a PERMANOVA approach

#parallel nodes (16 nodes).

# Results saved in 4 files at "/home/baldanzi/Sleep_apnea/Results/"
# File model 1 - "permanova_model1_t90_ad.tsv"
# File model 2 - "permanova_model2_t90_ad.tsv"
# File model 3 - "permanova_model3_t90_ad.tsv"
# File sensitivity analysis - permanova_model3_nomed_t90_ad.tsv

# Loading packages 
pacman::p_load(data.table, vegan, ggplot2,parallel)

rm(list = ls())
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.odi = fread('/home/baldanzi/Datasets/sleep_SCAPIS/validodi.MGS.Upp.tsv',
                  header=T, na.strings=c("", "NA")) 
setnames(valid.odi, "pob", "placebirth")
setnames(valid.odi, "sat90", "T90")
valid.odi = valid.odi[!is.na(T90),]

# Importing AD matrix 
print("Import AD matrix")
AD <-  fread('T90.ADmatrix.csv', header=T, sep = ',')
a <-  colnames(AD)
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
  result[1,"N"]=sum(a)
  result$variables = variables 
  result = result[,c("variables",names(result[,-ncol(result)]))]
  return(result)
}

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

# Making sure that BC and dades have the same order of observations 
dades = dades[match(rownames(AD),dades$SCAPISid),]

# Outcome - character name (length=1) with matrix distance 
outc = "AD"

# Main Exposure - character name (length=1)
expo = "T90"

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate + received 
model1 <-   c("age", "Sex", "Alkohol","smokestatus","plate")
# model 2 = model 1 + BMI 
model2 <-  c(model1,"BMI")
# model 3 = model 2 + fiber intake + physical activity + education + country of birth 
model3 <-  c(model2, "Fibrer","Energi_kcal", "leisurePA", "educat","placebirth")

# Runing PERMANOVA in parallel ####
set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outc","expo","dades","model1","model2","model3"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
res1 = PermanovaFunction(outcome = outc, exposure = expo, covari = model1, data = dades, nodes = nod)
res2 = PermanovaFunction(outcome = outc, exposure = expo, covari = model2, data = dades, nodes = nod)
res3 = PermanovaFunction(outcome = outc, exposure = expo, covari = model3, data = dades, nodes = nod)

stopCluster(cl)

# Saving result
fwrite(res1, file = paste0(output,"permanova_model1_t90_ad.tsv"), sep="\t")
fwrite(res2, file = paste0(output,"permanova_model2_t90_ad.tsv"), sep="\t")
fwrite(res3, file = paste0(output,"permanova_model3_t90_ad.tsv"), sep="\t")



#---------------------------------------------------------------------------#
# Sensitivity analysis excluding medication users 

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
fwrite(res3.nomed, file = paste0(output,"permanova_model3_nomed_t90_ad.tsv"), sep="\t")

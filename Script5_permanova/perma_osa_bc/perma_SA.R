# Sleep apnea and microbiota project 

# Gabriel Baldanzi - 2021-07-08

# Sensitivity analysis excluding medication users 

# Sensitivity analysis - remove individuals who use medication 
dades.sa <-  copy(valid.ahi)
dades.sa <-  dades.sa[ppi == "no",] #60
dades.sa <-  dades.sa[metformin == "no",] #57 
dades.sa <-  dades.sa[hypermed == "no",] #593
dades.sa <-  dades.sa[dyslipmed == "no",] # 239
nrow(dades.sa) #2318

# Transforming two-level factor variables into numeric variables 
a= c("Sex")
dades.sa[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades.sa[,a, with=F]))))]

# Transforming factor variables 
dades.sa[,plate:=as.factor(dades.sa$plate)]

# Making sure that BC and dades.sa have the same observations 
BC = as.data.frame(BC) # Transform BC from matrix to data.frame
cols = dades.sa[,SCAPISid] # Pass the dades.sa observations ids to a vector
BC = BC[,cols] # Only keep the BC columns that correspond to dades.sa observations
BC$SCAPISid = rownames(BC) # Creates a BC variable with the rownames(BC) 
BC = BC[BC$SCAPISid %in% dades.sa[,SCAPISid], ] # Exclude BC rows that are not present in dades.sa 
BC$SCAPISid = NULL # Exclude the SCAPISid column 

BC = as.matrix(BC) # Transform BC back to a matrix

dades.sa = dades.sa[match(rownames(BC),dades.sa$SCAPISid),]  # Makes that BC and dades.sa are in the same order

# Runing PERMANOVA in parallel ####
print("PERMANOVA OSA categoies and BC - Sensitivity Analysis")
print(" ")
set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outc","expo","dades.sa","SA"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
res = PermanovaFunction(outcome = outc, exposure = expo, covari = SA, data = dades.sa, nodes = nod)

stopCluster(cl)

# Saving results 
fwrite(res, file = paste0(output,"permanova_SA_osa_bc.tsv"), sep="\t")
rm(dades.sa)

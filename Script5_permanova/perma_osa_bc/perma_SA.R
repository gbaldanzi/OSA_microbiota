# Sleep apnea and microbiota project 

# Gabriel Baldanzi - 2021-07-08

# Sensitivity analysis excluding medication users 

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

# Making sure that BC and dades have the same observations 
BC = as.data.frame(BC) # Transform BC from matrix to data.frame
cols = dades[,SCAPISid] # Pass the dades observations ids to a vector
BC = BC[,cols] # Only keep the BC columns that correspond to dades observations
BC$SCAPISid = rownames(BC) # Creates a BC variable with the rownames(BC) 
BC = BC[BC$SCAPISid %in% dades[,SCAPISid], ] # Exclude BC rows that are not present in dades 
BC$SCAPISid = NULL # Exclude the SCAPISid column 

BC = as.matrix(BC) # Transform BC back to a matrix

dades = dades[match(rownames(BC),dades$SCAPISid),]  # Makes that BC and dades are in the same order

# Runing PERMANOVA in parallel ####
print("PERMANOVA OSA categoies and BC - Sensitivity Analysis")
print(" ")
set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outc","expo","dades","SA"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
res = PermanovaFunction(outcome = outc, exposure = expo, covari = SA, data = dades, nodes = nod)

stopCluster(cl)

# Saving results 
fwrite(res3.nomed, file = paste0(output,"permanova_SA_osa_bc.tsv"), sep="\t")

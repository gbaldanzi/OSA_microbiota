# Sleep apnea and microbiota project 

# Gabriel Baldanzi - 2021-07-08

# Sensitivity analysis excluding medication users 

# Sensitivity analysis - remove individuals who use medication 
dade.sa <-  copy(valid.t90)
dade.sa <-  dade.sa[ppi == "no",] #60
dade.sa <-  dade.sa[metformin == "no",] #57 
dade.sa <-  dade.sa[hypermed == "no",] #593
dade.sa <-  dade.sa[dyslipmed == "no",] # 239
nrow(dade.sa) #2576

# Transforming two-level factor variables into numeric variables 
a= c("Sex")
dade.sa[,(a):=as.data.frame(data.matrix(data.frame(unclass(dade.sa[,a, with=F]))))]


# Making sure that BC and dade.sa have the same observations 
  bcsa <-  as.data.frame(BC) # Transform BC from matrix to data.frame
  cols = dade.sa[,SCAPISid] # Pass the dade.sa observations ids to a vector
  bcsa = bcsa[,cols] # Only keep the bcsa columns that correspond to dade.sa observations
  bcsa$SCAPISid = rownames(bcsa) # Creates a bcsa variable with the rownames(bcsa) 
  bcsa = bcsa[bcsa$SCAPISid %in% dade.sa[,SCAPISid], ] # Exclude bcsa rows that are not present in dade.sa 
  bcsa$SCAPISid = NULL # Exclude the SCAPISid column 

bcsa = as.matrix(bcsa) # Transform bcsa back to a matrix

dade.sa = dade.sa[match(rownames(bcsa),dade.sa$SCAPISid),]  # Makes that bcsa and dade.sa are in the same order

# Outcome
outsa = "bcsa"

# Runing PERMANOVA in parallel ####
print("PERMANOVA T90 and bcsa - Sensitivity Analysis")
print(" ")
set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outsa","expo","dade.sa","SA","bcsa"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
res = PermanovaFunction(outcome = outsa, exposure = expo, covari = SA, data = dade.sa, nodes = nod)

stopCluster(cl)

# Saving results 
fwrite(res, file = paste0(output,"permanova_SA_t90_bc.tsv"), sep="\t")

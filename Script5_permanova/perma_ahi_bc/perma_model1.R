# Sleep apnea and microbiota project 

# Gabriel Baldanzi - 2021-07-08

set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outc","expo","dades","model1","model2","model3"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
print("PERMANOVA OSA categories and BC - Model1")
print(" ")
res = PermanovaFunction(outcome = outc, exposure = expo, covari = model1, data = dades, nodes = nod)
fwrite(res, file = paste0(output,"permanova_model1_osa_bc.tsv"), sep="\t")
stopCluster(cl)
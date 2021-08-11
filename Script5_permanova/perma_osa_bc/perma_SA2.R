# Sleep apnea and microbiota project 

# Gabriel Baldanzi - 2021-08-09

set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outc","expo","dades","model1","model2","model3","SA2"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
print("PERMANOVA OSA categories and BC - SA2")
print(" ")
res = PermanovaFunction(outcome = outc, exposure = expo, covari = SA2, data = dades, nodes = nod)
fwrite(res, file = paste0(output,"permanova_SA2_osa_bc.tsv"), sep="\t")
stopCluster(cl)
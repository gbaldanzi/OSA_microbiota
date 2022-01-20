# Sleep apnea and microbiota project 

# Gabriel Baldanzi - 2021-07-08

set.seed(123)
nod=16   # Number of workers to be used 
cl = makeCluster(nod)
clusterExport(cl, varlist = c("outc","expo","dades","basic.model"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
message(paste("PERMANOVA", expo, "and", outc, "- Basic Model"))
message(" ")
res = PermanovaFunction(outcome = outc, exposure = expo, covari = basic.model, data = dades, nodes = nod)
stopCluster(cl)
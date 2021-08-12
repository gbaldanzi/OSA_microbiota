# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-12

message("Linear model - AHI and MGS clr transformed") 

# CLR transformation 

clr.count <- fread('/home/baldanzi/Datasets/MGS/clean/MGS_clr.transformed_4839_upp.tsv')

a = c("SCAPISid",mgs.ahi)
clr.count <- clr.count[,a,with=F]

  dades <- copy(valid.ahi)

  dades[,grep("____",names(dades)):=NULL]

dades <- merge(dades, clr.count, by="SCAPISid", all.x=T, all.y=F)


# Linear regression 

res.ahi.clr <-  lapply(mgs.ahi, lm.mgs, phenotype = "ahi", covariates = model2, data = dades)

res.ahi.clr <- do.call(rbind,res.ahi.clr)
res.ahi.clr$transformation="clr"
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-12

message("Linear model - T90 and MGS clr transformed") 

# CLR transformation 

  clr.count <- fread('/home/baldanzi/Datasets/MGS/clean/MGS_clr.transformed_4839_upp.tsv')

  a = c("SCAPISid",mgs.t90)
  clr.count <- clr.count[,a,with=F]

  dades <- copy(valid.t90)

  dades[,grep("____",names(dades)):=NULL]

  dades <- merge(dades, clr.count, by="SCAPISid", all.x=T, all.y=F)


# Linear regression 

res.t90.clr <-  lapply(mgs.t90, lm.mgs, phenotype = "t90", covariates = model2, data = dades)

res.t90.clr <- do.call(rbind,res.t90.clr)
res.t90.clr$transformation="clr"
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-12

message("Linear model - T90 and MGS log+1 transformed") 

# Log+1 transformation 
dades = copy(valid.t90)
dades[, (mgs.t90) := lapply(.SD, function(x){log(x+1)}), .SDcols = mgs.t90]

# Linear regression 

res.t90.log <-  lapply(mgs.t90, lm.mgs, phenotype = "t90", covariates = model2, data = dades)

res.t90.log <- do.call(rbind,res.t90.log)
res.t90.log$transformation="log+1"
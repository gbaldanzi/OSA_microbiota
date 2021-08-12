# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-12

message("Linear model - AHI and MGS log+1 transformed") 

# Log+1 transformation 
dades = copy(valid.ahi)
dades[, (mgs.ahi) := lapply(.SD, function(x){log(x+1)}), .SDcols = mgs.ahi]

# Linear regression 

res.ahi.log <-  lapply(mgs.ahi, lm.mgs, phenotype = "ahi", covariates = model2, data = dades)

res.ahi.log <- do.call(rbind,res.ahi.log)
res.ahi.log$transformation="log+1"
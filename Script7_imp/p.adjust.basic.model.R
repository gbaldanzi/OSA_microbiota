# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Multiple comparison adjustment after imputation of AHI and 
# correlation with gut species using the basic model STATA 

library(data.table)

results.folder <- "/home/baldanzi/Sleep_apnea/Results/"

res.cor.basic.model1 <- fread(paste0(results.folder,"cor_ahi_imput_mgs_1.tsv"))
res.cor.basic.model2 <- fread(paste0(results.folder,"cor_ahi_imput_mgs_2.tsv"))
res.cor.basic.model3 <- fread(paste0(results.folder,"cor_ahi_imput_mgs_3.tsv"))
res.cor.basic.model4 <- fread(paste0(results.folder,"cor_ahi_imput_mgs_4.tsv"))

res.cor.basic.model <- rbind(res.cor.basic.model1, res.cor.basic.model2,
                             res.cor.basic.model3, res.cor.basic.model4)

res.cor.basic.model[,q_value := p.adjust(p_value, method="BH")]

setcolorder(res.cor.basic.model, c("MGS","exposure","rho","p_value","q_value","N"))

fwrite(res.cor.basic.model, paste0(results.folder,"cor_ahi_imput_mgs.tsv"))


## Save the species associated with either ODI, T90 or the imputed AHI 

  res <- fread(paste0(results.folder,"cor_all.var_mgs.tsv"))

  res[q.value>=.001,q_value := round(q.value,3)] 
  
  mgs.odi.t90 <- res[exposure %in% c("odi","t90") & q.value<.05 , MGS]
  
    cutlast <- function(char,n){
    l <- nchar(char)
    a <- l-n+1
    return(substr(char,a,l))
    }
    
    mgs.odi.t90 <- cutlast(mgs.odi.t90,9)
    
    mgs.odi.t90 <- gsub("[.]","_", mgs.odi.t90)
  
  
  res.cor.basic.model[q_value<0.001, qvalue := round(q_value,3)]
  
  mgs.ahi = res.cor.basic.model[q_value<0.05, MGS]
  
  sig.mgs <- data.frame(sigmgs=unique(c(mgs.odi.t90, mgs.ahi)))
  
  require(foreign)
  write.dta(sig.mgs, paste0(results.folder,"sig_mgs_basicmodel.dta"))

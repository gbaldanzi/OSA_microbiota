# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-09-14

# This code will create a table with the results for the SA for the relevant MGS

  library(data.table)
  library(tidyverse)

  # Importing results 
  res.m2 <-  fread("/home/baldanzi/Sleep_apnea/Results/cor2_all.var_mgs.tsv")
  res.sa <-  fread("/home/baldanzi/Sleep_apnea/Results/corsa_all.var_mgs.tsv")
  
  mgs.fdr <- readRDS("/home/baldanzi/Sleep_apnea/Results/mgs.m2.rds")
  names(mgs.fdr) <- c("ahi", "t90", "odi")
  
  # Round function 
  round.fun <- function(dt) {
    dt[q.value>=0.001, q.value:=round(q.value, digits = 3)]
    dt[p.value>=0.001, p.value:=round(p.value, digits = 3)]
    dt[,rho:=cor.coefficient]
    dt[, rho:=round(rho, digits = 3)]
  }
  
  res.m2 <- round.fun(res.m2)
  res.m2[,model:="full_model"]
  res.sa <- round.fun(res.sa)
  
  # Keep only the significant MGS
  
  compare.model.fun <- function(model.res, ref.res, expo, sig.mgs){
    
      a <- c("MGS","rho", "p.value","q.value","model","N")
      model.temp <- model.res[MGS %in% sig.mgs[[expo]] & exposure == expo, a, with = F]
     
      ref.temp <- res.m2[MGS %in% sig.mgs[[expo]] & exposure == expo, a, with = F]
     

      comparison <- rbind(ref.temp, model.temp) %>% 
      pivot_wider(id_cols = MGS, 
                  names_from=model, 
                  values_from = c(rho, p.value ,q.value, N))                                                     
    setDT(comparison)
    return(comparison)
  }
  
  list.res <- NULL
  
  list.res <- lapply(c("ahi", "t90", "odi"), compare.model.fun, 
                       model.res=res.sa,  ref.res=res.m2, sig.mgs=mgs.fdr)
  names(list.res) <- c("ahi","t90","odi")

  

  saveRDS(list.res,file="/home/baldanzi/Sleep_apnea/Results/table.res_sa_mgs.fdr.rds")
  
  
  
  
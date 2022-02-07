# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Lastest update: 2022-02-02

# This code will produce a table from the correlations of AHI, T90, ODI, and
# BMI with gut microbiota species using the basic model

  # MGSs identified with basic model
  mgs.m1 <- readRDS('/home/baldanzi/Sleep_apnea/Results/mgs.m1.rds') 
  
  # Results from basic model
  # res.basic.model <- fread(paste0(output,"cor_all.var_mgs.tsv"))
  
  res.list <- list(AHI = res.basic.model[exposure=="ahi",],
                   T90 = res.basic.model[exposure=="t90",],
                   ODI = res.basic.model[exposure=="odi",])
  
  res.list = lapply(res.list,Clean.Correlation.Results)
  
  
  res.df = do.call(rbind,res.list)
  
  # Restricting results to mgs.m1
  res.df[res.df$MGS %in% mgs.m1,]

  
  # From long to wide
  res.df <- res.df %>% pivot_wider(id_cols = MGS, names_from=exposure, 
                                   values_from= c("rho","p.value","q.value","N"))
  
  # Save table 
  saveRDS(res.df,file="/home/baldanzi/Sleep_apnea/Results/table.res1.rds")
  

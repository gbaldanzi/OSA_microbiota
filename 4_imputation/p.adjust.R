# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Multiple comparison adjustment after imputation of AHI and 
# correlation with gut species using the basic model STATA 


  library(data.table)

  results.folder <- "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/"

  res_1 <- fread(paste0(results.folder,"cor_ahi_imput_mgs_1.tsv"))
  res_2 <- fread(paste0(results.folder,"cor_ahi_imput_mgs_2.tsv"))
  res_3 <- fread(paste0(results.folder,"cor_ahi_imput_mgs_3.tsv"))
  res_4 <- fread(paste0(results.folder,"cor_ahi_imput_mgs_4.tsv"))
  
  res.impute <- rbind(res_1,res_2,res_3,res_4)

  res.impute[,q_value := p.adjust(p_value, method="BH")]

  setcolorder(res.impute, c("MGS","exposure","rho","p_value","q_value","N"))

  fwrite(res.impute, paste0(results.folder,"cor_ahi_imput_mgs.tsv"))

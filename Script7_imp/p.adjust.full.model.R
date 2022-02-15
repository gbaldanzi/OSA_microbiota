# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Multiple comparison adjustment after imputation of AHI and 
# correlation with gut species using the full model STATA 

library(data.table)

results.folder <- "/home/baldanzi/Sleep_apnea/Results/"

res.cor.full.model <- fread(paste0(results.folder,"cor2_ahi_imput_mgs.tsv"))

res.cor.full.model[,q_value := p.adjust(p_value, method="BH")]

setcolorder(res.cor.full.model, c("MGS","exposure","rho","p_value","q_value","N"))

fwrite(res.cor.full.model, paste0(results.folder,"cor2_ahi_imput_mgs.tsv"))

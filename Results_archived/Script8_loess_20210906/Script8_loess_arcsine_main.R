# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-12

# Create LOESS curves for the scatter plots on MGS and sleep variables 

# MGS are arcsine transformed

rm(list = ls())
pacman::p_load(data.table,ggplot2,ggpubr)

# input and output folders 
input1 = "/home/baldanzi/Sleep_apnea/Results/"
input2 = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

#output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"
output.plot = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

#source(paste0(input2,'lm.mgs.function.R'))

  # Import full data 
  valid.ahi <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/validsleep_MGS.shannon_Upp.rds")
  valid.t90 <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/valid.t90_MGS.shannon_Upp.rds")

  # Import results 
  res.ahi <- fread(paste0(input1,"cor2_ahi_mgs.tsv"))
  res.t90 <- fread(paste0(input1,"cor2_t90_mgs.tsv"))
  res.bmi <- fread(paste0(input1,"cor2_BMI_mgs.tsv"))

  res.ahi[,q.value:=round(q.value,3)]
  res.t90[,q.value:=round(q.value,3)]
  res.bmi[,q.value:=round(q.value,3)]

  # Select the relevant MGSs 
  mgs.bmi <- res.bmi[q.value<.05,]
  mgs.ahi <- res.ahi[q.value<.05 & !MGS %in% mgs.bmi$MGS,MGS]
  mgs.t90 <- res.t90[q.value<.05 & !MGS %in% mgs.bmi$MGS,MGS]
  
  a=c(mgs.ahi,"ahi")
  plotdata.ahi <- copy(valid.ahi[,a,with=F])
  a=c(mgs.t90,"t90")
  plotdata.t90 <- copy(valid.t90[,a,with=F])
  
  plotdata.ahi[,get("mgs.ahi"):=lapply(.SD,function(x){asin(sqrt(x/100))}), .SDcols=mgs.ahi]
  plotdata.t90[,get("mgs.t90"):=lapply(.SD,function(x){asin(sqrt(x/100))}), .SDcols=mgs.t90]
  
  #LOESS plot 
  source("Script8_loess/Script8_loess_arcsine_ahi.R")
  source("Script8_loess/Script8_loess_arcsine_t90.R")
  
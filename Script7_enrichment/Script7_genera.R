# Script 7 - Enrichment analysis using MGS identified with model 2

# Gabriel Baldanzi v1 2021-06-22

# Loading packages 
pacman::p_load(data.table,ggplot2, tidyr, fgsea)

source('MGS.Enrich.function.R')

# input and output folders 
input = "/home/baldanzi/Sleep_apnea/Results/"
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# List of genera ####
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  
  gg = unique(taxonomy[,genus])
  
  list.genera = lapply(gg,function(i){
    taxonomy[genus==i,maintax_mgs]
  })
  
  names(list.genera)=gg


# Importing results 
  res.ahi <- fread(paste0(input,"cor_ahi_mgs.tsv"))
  res.bmi <- fread(paste0(input,"cor_BMI_mgs.tsv"))
  res.t90 <- fread(paste0(input,"cor_t90_mgs.tsv"))

  list.res = list(res.ahi,res.bmi,res.t90)

  # Enrichment analysis 
  
  res.both <-  lapply(list.res,MGS.Enrich.Analysis,  
                     p.value.name="p.value",
                     cor.var.name = "cor.coefficient",
                     MGS.var.name = "MGS",
                     enrich.var.list = list.genera,
                     direction = "both")
  names(res.both) <- c("ahi","bmi","t90")
  for(i in 1:3){setnames(res.both[[i]],"pathway","genera")}
  saveRDS(res.both, file = paste0(output,"ea_genera_both.rds"))

  res.pos <-  lapply(list.res,MGS.Enrich.Analysis,  
                    p.value.name="p.value",
                    cor.var.name = "cor.coefficient",
                    MGS.var.name = "MGS",
                    enrich.var.list = list.genera,
                    direction = "positive")
  names(res.pos) <- c("ahi","bmi","t90")
  for(i in 1:3){setnames(res.pos[[i]],"pathway","genera")}
  saveRDS(res.pos, file = paste0(output,"ea_genera_pos.rds"))
  
  res.neg  <- lapply(list.res,MGS.Enrich.Analysis,  
                     p.value.name="p.value",
                     cor.var.name = "cor.coefficient",
                     MGS.var.name = "MGS",
                     enrich.var.list = list.genera,
                     direction = "negative")
  names(res.neg) <- c("ahi","bmi","t90")
  for(i in 1:3){setnames(res.neg[[i]],"pathway","genera")}
  saveRDS(res.neg, file = paste0(output,"ea_genera_neg.rds"))
  


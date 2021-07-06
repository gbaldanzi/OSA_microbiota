# Script 7 - Enrichment analysis using MGS identified with model 1

#Rho - ea based on rho correlation coefficient instead of p.value

# Gabriel Baldanzi 2021-07-05

source('plot.ea.top.function.R')

# Loading packages 
pacman::p_load(data.table,ggplot2, tidyr, fgsea)

# input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results
  res.ea <- readRDS(file = paste0(input,"ea_genera_rho.rds"))
  
  res.ahi <- fread(paste0(input,"cor_ahi_mgs.tsv"))
  res.bmi <- fread(paste0(input,"cor_BMI_mgs.tsv"))
  res.t90 <- fread(paste0(input,"cor_t90_mgs.tsv"))
  
  list.res = list(res.ahi,res.bmi,res.t90)
  
  # Importing data 

# Genera list 
# List of genera ####
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")
  
  gg = unique(taxonomy[,genus])

  list.genera = lapply(gg,function(i){ taxonomy[ genus == i , maintax_mgs ] } )

  names(list.genera)=gg

# Top genera plot 
  g.ahi=    plot.ea.top(data = list.res[[1]],
                ea.res = res.ea[[1]],
                pathway.name="genera",
                cor.var.name = "cor.coefficient",
                MGS.var.name = "MGS",
                enrich.var.list = list.genera,
                output.plot=output.plot)
  
  save(g.ahi, file = paste0(output.plot,"ea_topplot_ahi"))
      
  g.bmi=    plot.ea.top(data = list.res[[2]],
                        ea.res = res.ea[[2]],
                        pathway.name="genera",
                        cor.var.name = "cor.coefficient",
                        MGS.var.name = "MGS",
                        enrich.var.list = list.genera,
                        output.plot=output.plot)
  
  save(g.bmi, file = paste0(output.plot,"ea_topplot_bmi"))
  
  g.t90=    plot.ea.top(data = list.res[[3]],
                        ea.res = res.ea[[3]],
                        pathway.name="genera",
                        cor.var.name = "cor.coefficient",
                        MGS.var.name = "MGS",
                        enrich.var.list = list.genera,
                        output.plot=output.plot)
  
  save(g.t90, file = paste0(output.plot,"ea_topplot_t90"))
    



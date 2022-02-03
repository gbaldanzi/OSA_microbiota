# Script 9 - Enrichment analysis using basic model correlation coefficients 

# Gabriel Baldanzi 

# Last update: 2022-02-03

rm(list=ls())

# Loading packages 
pacman::p_load(data.table,ggplot2, tidyr, fgsea)

  # Functions 
  source("Script0_functions/MGS.Enrich.function.R")

# input and output folders 
  results.folder = "/home/baldanzi/Sleep_apnea/Results/"
  input = "/home/baldanzi/Datasets/MGS/original/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Import results from basic model 
  res <- fread(paste0(results.folder,"cor_all.var_mgs.tsv"))
  
  res[,mgs:=cutlast(MGS,9)]
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   ODI = res[exposure=="odi",],
                   BMI = res[exposure=="BMI",])
  

  # List of modules ####
  load(paste0(input,'MGS_HG3A.GMMs2MGS.RData')) # object = MGS_HG3A.GMMs2MGS
  list.modules <-  MGS_HG3A.GMMs2MGS

  # Enrichment analysis   
    # Positive correlations ####
    message("Positive correlations")
    res.pos = list()[1:4]
    
    for(i in 1:4){
    
      res.pos[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                       p.value.name="p.value",
                       cor.var.name = "rho",
                       MGS.var.name = "mgs",
                       enrich.var.list = list.modules,
                       direction = "positive", 
                       maxSize = Inf)
          names(res.pos)[i] <- names(res.list)[i]
    }

    
      # Save results from enrichment analysis in the positive correlations  
    
        res.pos <- do.call(rbind,res.pos)
        fwrite(res.pos, file = paste0(results.folder,"ea_GMM_pos.tsv"))
      
      # Negative correlations ####

        message("Negative correlations")
        # Removing the single MGS that is negative correlated and not present in list of modules
       # res.list <- lapply(res.list,function(x){x[mgs!="HG3A.1213",]})
        
        res.neg = list()[1:4]
        
        for(i in 1:4){
          
        res.neg[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                           p.value.name="p.value",
                           cor.var.name = "rho",
                           MGS.var.name = "mgs",
                           enrich.var.list = list.modules,
                           direction = "negative", 
                           maxSize = Inf)
        names(res.neg)[i] <- names(res.list)[i]
        }
        
        # Save results from enrichment analysis with negative correlations 
        res.neg <- do.call(rbind, res.neg)
        fwrite(res.neg, file = paste0(results.folder,"ea_GMM_neg.tsv"))
        
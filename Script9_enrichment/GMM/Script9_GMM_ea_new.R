# Script 9 - Enrichment analysis using main model correlation coefficients 

# Gabriel Baldanzi 

# Last update: 2022-02-24

rm(list=ls())

# Loading packages 
pacman::p_load(data.table,ggplot2, tidyr, fgsea,stringr)

  # Functions 
  source("Script0_functions/MGS.Enrich.function.R")

# input and output folders 
  results.folder = "/home/baldanzi/Sleep_apnea/Results/"
  input = "/home/baldanzi/Datasets/MGS/original/"
  
  # Import pathways/modules names 
  gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  gmm.names[,Name:=str_to_title(Name)]
  gmm.names[,Name:=gsub("Ii","II",Name)]
  
  # Import results from basic model 
  res <- fread(paste0(results.folder,"cor.bmi_all.var_mgs.tsv"))
  
  res[,mgs:=cutlast(MGS,9)]
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   ODI = res[exposure=="odi",])
  

  # List of modules ####
  load(paste0(input,'MGS_HG3A.GMMs2MGS.RData')) # object = MGS_HG3A.GMMs2MGS
  list.modules <-  MGS_HG3A.GMMs2MGS
  
  for(m in names(list.modules)){
    list.modules[[m]] <- list.modules[[m]][list.modules[[m]] %in% res[rho>0,mgs]]
  }

  # Enrichment analysis   
    # Positive correlations ####
    message("Positive correlations")
    res.pos = list()[1:3]
    
    for(i in 1:3){
    
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
        res.pos <- merge(res.pos,gmm.names,by.x="pathway",by.y="Module",all.x=T,all.y=F)
        fwrite(res.pos, file = paste0(results.folder,"ea_GMM_pos.tsv"))
      
      # Negative correlations ####

        message("Negative correlations")

        list.modules <-  MGS_HG3A.GMMs2MGS
        
        for(m in names(list.modules)){
          list.modules[[m]] <- list.modules[[m]][list.modules[[m]] %in% res[rho<0,mgs]]
        }
        
        res.neg = list()[1:3]
        
        for(i in 1:3){
          
        res.neg[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                           p.value.name="p.value",
                           cor.var.name = "rho",
                           MGS.var.name = "mgs",
                           enrich.var.list = list.modules,
                           direction = "negative", 
                           maxSize = 700)
        names(res.neg)[i] <- names(res.list)[i]
        }
        
        # Save results from enrichment analysis with negative correlations 
        res.neg <- do.call(rbind, res.neg)
        res.neg <- merge(res.neg,gmm.names,by.x="pathway",by.y="Module",all.x=T,all.y=F)
        fwrite(res.neg, file = paste0(results.folder,"ea_GMM_neg.tsv"))
        
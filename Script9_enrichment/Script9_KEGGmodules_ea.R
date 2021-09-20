# Script 7 - Enrichment analysis using model2 correlation coefficients 

# Gabriel Baldanzi v1 2021-06-30

# Loading packages 
pacman::p_load(data.table,ggplot2, tidyr, fgsea,rio)

source("MGS.Enrich.function.R")

# input and output folders 
  input1 = "/home/baldanzi/Sleep_apnea/Results/"
  input2 = "/home/baldanzi/Datasets/MGS/original/"
  output = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Import results 
  res <- fread(paste0(input1,"cor_all.var_mgs.tsv"))
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   BMI = res[exposure=="BMI",])

# List of modules ####
  load(paste0(input2,'MGS_HG3A.keggModule2MGS.RData')) # object = MGS_HG3A.keggModule2MGS
  list.modules = MGS_HG3A.keggModule2MGS
  
  modules=import(paste0(input2,"upugut03.keggModuleComp.percent.tsv")) # modules number and description 
  modules <- modules[,1:4]

  # Enrichment analysis   
    # Positive correlations 
    message("Positive correlations")
    res.pos = list()[1:3]
    
    for(i in 3:1){
    print(i)
      
      res.pos[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                       p.value.name="p.value",
                       cor.var.name = "cor.coefficient",
                       MGS.var.name = "mgs",
                       enrich.var.list = list.modules,
                       direction = "positive")
          names(res.pos)[i] <- names(res.list)[i]
    }

    
      # Merging with module annotation 
      
        res.pos <- lapply(res.pos,function(x){merge(x,modules, by.x="pathway",by.y="Module")})
    
        res.pos <- do.call(rbind,res.pos)
        fwrite(res.pos, file = paste0(output,"ea_modules_pos.tsv"))
      
    # Negative correlations

        message("Negative correlations")
        # Removing the single MGS that is negative correlated and not present in list of modules
        res.list <- lapply(res.list,function(x){x[mgs!="HG3A.1213",]})
        res.neg = list()[1:3]
        
        for(i in 1:3){
          print(i)
        res.neg[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                           p.value.name="p.value",
                           cor.var.name = "cor.coefficient",
                           MGS.var.name = "mgs",
                           enrich.var.list = list.modules,
                           direction = "negative")
        names(res.neg)[i] <- names(res.list)[i]
        }
        
        res.neg <- lapply(res.neg,function(x){merge(x,modules, by.x="pathway",by.y="Module")})
        
        res.neg <- do.call(rbind, res.neg)
        fwrite(res.neg, file = paste0(output,"ea_modules_neg.tsv"))
        
    # Both correlations 
        message("Both directions")
        res.both = list()[1:3]
        for(i in 1:3){
          print(i)
        res.both[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                           p.value.name="p.value",
                           cor.var.name = "cor.coefficient",
                           MGS.var.name = "mgs",
                           enrich.var.list = list.modules,
                           direction = "both")
        names(res.both)[i] <- names(res.list)[i]
        }
        
        # Merging with module annotation 
        res.both <- lapply(res.both,function(x){merge(x,modules, by.x="pathway",by.y="Module")})
        
        res.both <- do.call(rbind,res.both)
        fwrite(res.both, file = paste0(output,"ea_modules_both.tsv"))

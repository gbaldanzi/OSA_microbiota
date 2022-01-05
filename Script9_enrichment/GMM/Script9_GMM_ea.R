# Script 9 - Enrichment analysis using model2 correlation coefficients 

# Gabriel Baldanzi v1 2021-10-25

# Last update: 

# Loading packages 
pacman::p_load(data.table,ggplot2, tidyr, fgsea,rio, ggvenn)

source("MGS.Enrich.function.R")

# input and output folders 
  input1 = "/home/baldanzi/Sleep_apnea/Results/"
  input2 = "/home/baldanzi/Datasets/MGS/original/"
  output = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Import results 
  res <- fread(paste0(input1,"cor_all.var_mgs.tsv"))
  
  
  cutlast <- function(char,n){
    l <- nchar(char)
    a <- l-n+1
    return(substr(char,a,l))
  }
  
  res[,mgs:=cutlast(MGS,9)]
  
  
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   ODI = res[exposure=="odi",],
                   BMI = res[exposure=="BMI",])

# List of modules ####
  load(paste0(input2,'MGS_HG3A.GMMs2MGS.RData')) # object = MGS_HG3A.GMMs2MGS
  list.modules <-  MGS_HG3A.GMMs2MGS

  # Enrichment analysis   
    # Positive correlations 
    message("Positive correlations")
    res.pos = list()[1:4]
    
    for(i in 1:4){
    print(i)
      
      res.pos[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                       p.value.name="p.value",
                       cor.var.name = "cor.coefficient",
                       MGS.var.name = "mgs",
                       enrich.var.list = list.modules,
                       direction = "positive", 
                       maxSize = Inf)
          names(res.pos)[i] <- names(res.list)[i]
    }

    
      # Merging with module annotation 
    
        res.pos <- do.call(rbind,res.pos)
        fwrite(res.pos, file = paste0(output,"ea_GMM_pos.tsv"))
      
    # Negative correlations

        message("Negative correlations")
        # Removing the single MGS that is negative correlated and not present in list of modules
       # res.list <- lapply(res.list,function(x){x[mgs!="HG3A.1213",]})
        res.neg = list()[1:4]
        
        for(i in 1:4){
          print(i)
        res.neg[[i]] <-  MGS.Enrich.Analysis(res.list[[i]],  
                           p.value.name="p.value",
                           cor.var.name = "cor.coefficient",
                           MGS.var.name = "mgs",
                           enrich.var.list = list.modules,
                           direction = "negative", 
                           maxSize = Inf)
        names(res.neg)[i] <- names(res.list)[i]
        }
        
        
        res.neg <- do.call(rbind, res.neg)
        fwrite(res.neg, file = paste0(output,"ea_GMM_neg.tsv"))
        
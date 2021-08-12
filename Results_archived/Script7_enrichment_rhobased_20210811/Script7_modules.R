# Script 7 - Enrichment analysis using MGS identified with model 1

# Gabriel Baldanzi v1 2021-06-30

# Loading packages 
pacman::p_load(data.table,ggplot2, tidyr, fgsea,rio)

source("MGS.Enrich.function.R")

# input and output folders 
input1 = "/home/baldanzi/Sleep_apnea/Results/"
input2 = "/home/baldanzi/Datasets/MGS/original/"
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing results 
res.ahi <- fread(paste0(input1,"cor_ahi_mgs.tsv"))
res.bmi <- fread(paste0(input1,"cor_BMI_mgs.tsv"))
res.t90 <- fread(paste0(input1,"cor_t90_mgs.tsv"))


list.res = list(res.bmi,res.ahi,res.t90)
names(list.res) = c("bmi","ahi","t90")


# List of modules ####
load(paste0(input2,'MGS_HG3A.keggModule2MGS.RData')) # object = MGS_HG3A.keggModule2MGS
list.modules = MGS_HG3A.keggModule2MGS
modules=import(paste0(input2,"upugut03.keggModuleComp.percent.tsv")) # modules number and description 
    modules=modules[,1:4]

  # Enrichment analysis   
    # Positive correlations 
    print("Positive correlations")
    res.pos = list()[1:3]
    
    #for(i in 1:3){
     # print(i)
    print(1)
      res.pos[[1]] <-  MGS.Enrich.Analysis(list.res[[1]],  
                       p.value.name="p.value",
                       cor.var.name = "cor.coefficient",
                       MGS.var.name = "mgs",
                       enrich.var.list = list.modules,
                       direction = "positive")
      print(2)
      res.pos[[2]] <-  MGS.Enrich.Analysis(list.res[[2]],  
                                           p.value.name="p.value",
                                           cor.var.name = "cor.coefficient",
                                           MGS.var.name = "mgs",
                                           enrich.var.list = list.modules,
                                           direction = "positive")
      print(3)
      res.pos[[3]] <-  MGS.Enrich.Analysis(list.res[[3]],  
                                           p.value.name="p.value",
                                           cor.var.name = "cor.coefficient",
                                           MGS.var.name = "mgs",
                                           enrich.var.list = list.modules,
                                           direction = "positive")
    
      names(res.pos) <- c("bmi","ahi","t90")
    
      # Merging with module annotation 
      
        res.pos <- lapply(res.pos,function(x){merge(x,modules, by.x="pathway",by.y="Module")})
    
        saveRDS(res.pos, file = paste0(output,"ea_modules_pos.rds"))
      
    # Negative correlations

        print("Negative correlations")
        # Removing the single MGS that is negative correlated and not present in list of modules
        list.res <- lapply(list.res,function(x){x[mgs!="HG3A.1213",]})
        res.neg = list()[1:3]
        
        for(i in 1:3){
          print(i)
        res.neg[[i]] <-  MGS.Enrich.Analysis(list.res[[i]],  
                           p.value.name="p.value",
                           cor.var.name = "cor.coefficient",
                           MGS.var.name = "mgs",
                           enrich.var.list = list.modules,
                           direction = "negative")
        }
        names(res.neg) <- c("bmi","ahi","t90")
        res.neg <- lapply(res.pos,function(x){merge(x,modules, by.x="pathway",by.y="Module")})
        
        saveRDS(res.neg, file = paste0(output,"ea_modules_neg.rds"))
        
    # Both correlations 
        print("Both directions")
        res.both = list()[1:3]
        for(i in 1:3){
          print(i)
        res.both[[i]] <-  MGS.Enrich.Analysis(list.res[[i]],  
                           p.value.name="p.value",
                           cor.var.name = "cor.coefficient",
                           MGS.var.name = "mgs",
                           enrich.var.list = list.modules,
                           direction = "both")
        }
        names(res.both) <- c("bmi","ahi","t90")
        
        # Merging with module annotation 
        
        res.both <- lapply(res.both,function(x){merge(x,modules, by.x="pathway",by.y="Module")})
        
        saveRDS(res.both, file = paste0(output,"ea_modules_both.rds"))

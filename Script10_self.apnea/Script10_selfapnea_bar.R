  
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-09-16

# Last update: 2021-09-16

# To plot prevalence of MGS identified with fully adjusted model according to 
#self reported sleep apnea 

  library(data.table) 
  library(tidyverse)
  library(vegan)

  pheno=fread("/home/baldanzi/Datasets/sleep_SCAPIS/SCAPIS-DATA-PETITION-170-20210315.csv", header = T, na.strings=c("", "NA"))
  setnames(pheno, "Subject", "SCAPISid") #renaming the id variable

  # Uploading MGS data ####
  data.MGS = fread("/home/baldanzi/Datasets/MGS/clean/MGS_relative_abundance_4839_upp_4980_malmo.tsv")
  data.MGS[SCAPISid=="", SCAPISid:=id]

  pheno <- merge(pheno, data.MGS, by="SCAPISid", all.y=T)
  
  # Recoding the self reported sleep apnea variable 
  pheno[cqhe058=="NO_ANSWER", cqhe058:=NA]
  
  # Identified MGS in fully-adjusted model 
  res.m2 <-  fread("/home/baldanzi/Sleep_apnea/Results/cor2_all.var_mgs.tsv")
    res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]
    mgs.bmi <- res.m2[q.value<0.05 & exposure=="BMI",MGS]
    mgs.fdr <- unique(res.m2[q.value<0.05 & exposure %in% c("ahi","t90"),MGS])
    mgs.fdr <- mgs.fdr[!mgs.fdr %in% mgs.bmi] # The identified MGS #
  
    
    # Determine direction of correlation between MGS and sleep apnea 
    cor_mgs.fdr <- res.m2[MGS %in% mgs.fdr & exposure %in% c("ahi","t90"),
                        .(MGS,exposure,cor.coefficient,p.value,q.value)]
  
    cor_mgs.fdr <- cor_mgs.fdr %>% pivot_wider(id_cols=MGS, names_from = exposure, 
                                             values_from = c(cor.coefficient,p.value,q.value))
  
    setDT(cor_mgs.fdr)
  
    all(cor_mgs.fdr[cor.coefficient_t90>0,cor.coefficient_ahi]>0) #TRUE
    all(cor_mgs.fdr[cor.coefficient_t90<0,cor.coefficient_ahi]<0) #TRUE
  
    cor_mgs.fdr[cor.coefficient_t90>0,direction:="pos"]
    cor_mgs.fdr[cor.coefficient_t90<0,direction:="neg"]
    
    mgs.pos <- cor_mgs.fdr[direction=="pos",MGS]
    mgs.neg <- cor_mgs.fdr[direction=="neg",MGS]
    
    # Presence absence 
    pa <- paste0("pa_",mgs.fdr)
    pheno[,(pa):= lapply(.SD, decostand, method="pa"), .SDcols = mgs.fdr]
    
    # Plot function 
    bar.plot.prevalence.fun <- function(y.axis,group,data){
      require(data.table)
      require(ggplot2)
      require(dplyr)
      dades <- data[!is.na(get(group)),]
      
      setDF(dades)
      
      dades <- dades %>% group_by(get(group)) %>% summarise_at(pa,function(x){sum(x)/length(x)})
      names(dades)[1] <- group
      names(dades) <- gsub("pa_","",names(dades))
      
      ggplot(dades) + geom_bar(aes_string(x=group, y=y.axis), stat = 'identity') + 
        ggtitle(y.axis) + 
        xlab("Self-reported sleep apnea") + ylab(y.axis) +
        theme_classic() +
        theme(plot.title = element_text(hjust = .5, face = 'bold'),
              plot.subtitle = element_text(hjust = .5, face = 'bold', size=16),
              axis.title.x = element_text(size=14),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 10)) 
      
    }
    
    # Creating the plots 
    if(dir.exists("Self_apnea_mgs/pos.mgs")==FALSE){dir.create("Self_apnea_mgs/pos.mgs")}
    if(dir.exists("Self_apnea_mgs/neg.mgs")==FALSE){dir.create("Self_apnea_mgs/neg.mgs")}
    
    
    plots.pos <- lapply(mgs.pos, 
                        bar.plot.prevalence.fun, 
                        group = "cqhe058",
                        data = pheno)
    
    names(plots.pos) <- mgs.pos
    
    plots.neg <- lapply(mgs.neg, 
                        bar.plot.prevalence.fun, 
                        group = "cqhe058",
                        data = pheno)
    
    names(plots.neg) <- mgs.neg
    
    # Saving plots 
    
    setwd("/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/")
    
    
    for(i in 1:length(plots.pos)){
      ggsave(file = paste0("selfapnea_",names(plots.pos[i]),".png"), plot=plots.pos[[i]],
             path = './Self_apnea_mgs/pos.mgs/', dpi = 60)
    }
    
    for(i in 1:length(plots.neg)){
      ggsave(file = paste0("selfapnea_",names(plots.neg[i]),".png"), plot=plots.neg[[i]],
             path = './Self_apnea_mgs/neg.mgs/', dpi = 60)
    }
    
    
    
    
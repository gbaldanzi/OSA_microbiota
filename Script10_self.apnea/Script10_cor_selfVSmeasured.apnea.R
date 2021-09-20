# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Version 1: 2021-09-17
# Last Update: 2021-09-17


  rm(list = ls())

  library(data.table)
  library(tidyverse)
  
  # Import fully adjusted model results 
  res.m2 <-  fread("/home/baldanzi/Sleep_apnea/Results/cor2_all.var_mgs.tsv")
  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  mgs.bmi <- res.m2[q.value<0.05 & exposure=="BMI",MGS]
  res.m2.ahi <- res.m2[q.value<0.05 & exposure=="ahi" & !MGS %in% mgs.bmi, ]
  res.m2.t90 <- res.m2[q.value<0.05 & exposure=="t90" & !MGS %in% mgs.bmi, ]
  
  mgs.fdr <- unique(res.m2[q.value<0.05 & exposure %in% c("ahi","t90"),MGS])
  mgs.fdr <- mgs.fdr[!mgs.fdr %in% mgs.bmi] # The identified MGS #
  
  # Determine direction of correlation between MGS and sleep apnea 
  cor_mgs.fdr <- res.m2[MGS %in% mgs.fdr & exposure %in% c("ahi","t90"),
                        .(MGS,exposure,cor.coefficient,p.value,q.value)]
  
  cor_mgs.fdr <- cor_mgs.fdr %>% pivot_wider(id_cols=MGS, names_from = exposure, 
                                             values_from = c(cor.coefficient,p.value,q.value))
  
  mgs.pos <- cor_mgs.fdr[cor_mgs.fdr$cor.coefficient_t90>0,"MGS"]
  mgs.neg <- cor_mgs.fdr[cor_mgs.fdr$cor.coefficient_t90<0,"MGS"]
  
  # Import Self.Apnea correlations 
  res.self <- fread("/home/baldanzi/Sleep_apnea/Results/cor_self.apnea_mgs.tsv")
  
  n.self <- unique(res.self$N)
  
  res.self.ahi <- res.self[MGS %in% res.m2.ahi$MGS,]
  cor.ahi <- cor.test(res.self.ahi$cor.coefficient, res.m2.ahi$cor.coefficient, method="spearman")
  
  
  res.self.t90 <- res.self[MGS %in% res.m2.t90$MGS,]
  cor.t90 <- cor.test(res.self.t90$cor.coefficient, res.m2.t90$cor.coefficient, method="spearman")
  
  n.self <- unique(res.self$N)
  n.ahi <- unique(res.m2.ahi$N)
  n.t90 <- unique(res.m2.t90$N)
  
  # Scatter plots 
  p1 <- ggplot(data.frame(self.reported=res.self.ahi$cor.coefficient, 
                          measured.apnea=res.m2.ahi$cor.coefficient),
               aes(x=measured.apnea, y = self.reported)) +
    geom_point() + geom_smooth(method='lm',se=F) +
    ggtitle("MGS correlations for Self reported VS measured AHI", 
           subtitle = paste0("rho=",round(cor.ahi$estimate,2),
                             ", p=",round(cor.ahi$p.value,4))) +
           xlab(paste0("Cor. Coefficient - AHI and MGS\nn=",n.ahi)) +
           ylab(paste0("Cor. Coefficient - reported apnea and MGS\nn=",n.self))
  
  ggsave(file="scatter_self_ahi.png", plot=p1, dpi = 150, 
         path='/home/baldanzi/Sleep_apnea/Results/Plots/')
  
  
  p2 <- ggplot(data.frame(self.reported=res.self.t90$cor.coefficient, 
                          measured.apnea=res.m2.t90$cor.coefficient),
               aes(x=measured.apnea, y = self.reported))+
    geom_point() + geom_smooth(method='lm',se=F) +
    ggtitle("MGS correlations for Self reported VS measured T90", 
            subtitle = paste0("rho=",round(cor.t90$estimate,2),
                              ", p=",round(cor.t90$p.value,4))) +
            xlab(paste0("Cor. Coefficient - T90 and MGS\nn=",n.t90)) +
    ylab(paste0("Cor. Coefficient - reported apnea and MGS\nn=",n.self))
  
  ggsave(file="scatter_self_t90.png", plot=p2, dpi = 150, 
         path='/home/baldanzi/Sleep_apnea/Results/Plots/')
  
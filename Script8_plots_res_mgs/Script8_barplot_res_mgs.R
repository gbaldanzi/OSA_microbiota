# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-09-06

# Last update: 2021-09-20

# To investigate non-linear associations between MGS and AHI


rm(list = ls())
pacman::p_load(data.table,tidyverse,Hmisc, vegan)

# input and output folders 
input1 = "/home/baldanzi/Sleep_apnea/Results/"
input2 = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

#output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"
output.plot = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

  # Import full data 
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

  # Import results 
  res.m2 <- fread(paste0(input1,"cor2_all.var_mgs.tsv"))
  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]

  # Select the relevant MGSs 
  mgs.bmi <- res.m2[exposure =="BMI" & q.value<.05,MGS] # MGS correlated to BMI
  res.ahi <- res.m2[exposure =="ahi" & q.value<.05,] %>% filter(!MGS %in% mgs.bmi)
  mgs.ahi <- res.ahi$MGS
  res.t90 <- res.m2[exposure =="t90" & q.value<.05,] %>% filter(!MGS %in% mgs.bmi)
  mgs.t90 <- res.t90$MGS

  # Result list 
  res.list <- list(ahi = res.ahi, t90 = res.t90)

  res.list <- lapply(res.list, function(x){x[,rho:=paste0("rho = ",round(cor.coefficient,3))]})
  
  #Categories
  pheno[valid.t90=="yes" & t90==0, t90cat := 0 ]
  pheno[valid.t90=="yes" & t90!=0, T90cat :=  cut(t90,breaks = quantile(t90,
                                                                  probs = seq(0,1,by=.33),
                                                                  na.rm=T), include.lowest = T)]
  
  
  # Bar plots of prevalence 
  noms=unique(c(mgs.ahi, mgs.t90))
  pa <- paste0("pa_",noms)
  pheno[,(pa):= lapply(.SD, decostand, method="pa"), .SDcols = noms]
  
 
  # Function for plots 
  bar.plot.prevalence.fun <- function(y.axis,exposure,group,data){
    require(data.table)
    require(ggplot2)
    require(dplyr)
    require(Hmisc)
    dades <- data[!is.na(get(group)),]
  
    setDF(dades)
  
    dades <- dades %>% group_by(get(group)) %>% summarise_at(pa,function(x){sum(x)/length(x)})
    names(dades)[1] <- group
    names(dades) <- gsub("pa_","",names(dades))
  
    r <- res.list[[exposure]][MGS==y.axis,rho]
  
  
  
    ggplot(dades) + geom_bar(aes_string(x=group, y=y.axis), stat = 'identity') + 
      ggtitle(y.axis, subtitle = r) + 
      xlab(group) + ylab(y.axis) +
      theme_classic() +
      theme(plot.title = element_text(hjust = .5, face = 'bold'),
          plot.subtitle = element_text(hjust = .5, face = 'bold', size=16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 10)) 

  }

  # Creating the plots 
  dir.create("bar_ahi")
  
  
  message("Plot OSAcat")
  plots.ahi <- lapply(mgs.ahi, 
                    bar.plot.prevalence.fun, 
                    exposure = "ahi",
                    group = "OSAcat",
                    data = pheno[valid.ahi=='yes',])

  names(plots.ahi) <- mgs.ahi
  saveRDS(plots.ahi, file = paste0(output.plot,'/bar_ahi/barplot_ahi.rds'))

  for(i in 1:length(plots.ahi)){
      ggsave(file = paste0("bar_ahi_",names(plots.ahi[i]),".png"), plot=plots.ahi[[i]],
         path = './bar_ahi/')
      }

  dir.create("bar_t90")
  
  plots.t90 <- lapply(mgs.t90, 
                    bar.plot.prevalence.fun, 
                    group= "t90cat",
                    exposure = "t90",
                    data = pheno[valid.t90=='yes',])

  names(plots.t90) <- mgs.t90
  saveRDS(plots.t90, file = paste0(output.plot,'/bar_t90/barplot_t90.rds'))

  for(i in 1:length(plots.t90)){
  
    ggsave(file = paste0("bar_t90_",names(plots.t90[i]),".png"), plot=plots.t90[[i]],
         path = './bar_t90/')
  }


# plots of quantiles 
  plots.ahi <- lapply(mgs.ahi, 
                    bar.plot.prevalence.fun, 
                    group = "AHIquantile",
                    exposure = "ahi",
                    data = pheno[valid.ahi=='yes',])

  names(plots.ahi) <- mgs.ahi
  
  saveRDS(plots.ahi, file = paste0(output.plot,'/bar_ahi/barplot_quant_ahi.rds'))

  for(i in 1:length(plots.ahi)){
  
  ggsave(file = paste0("bar_quant_ahi_",names(plots.ahi[i]),".png"), plot=plots.ahi[[i]],
         path = './bar_ahi/')
  }


  plots.t90 <- lapply(mgs.t90, 
                    bar.plot.prevalence.fun, 
                    exposure = "t90", 
                    group= "T90quantile",
                    data = pheno[valid.t90=='yes',])

  names(plots.t90) <- mgs.t90
  
  saveRDS(plots.t90, file = paste0(output.plot,'/bar_t90/barplot_quant_t90.rds'))
  
  for(i in 1:length(plots.t90)){
  
      ggsave(file = paste0("bar_quant_t90_",names(plots.t90[i]),".png"), plot=plots.t90[[i]],
         path = './bar_t90/')
  }

  source('Script8_minibarplots.R')
  
  
  
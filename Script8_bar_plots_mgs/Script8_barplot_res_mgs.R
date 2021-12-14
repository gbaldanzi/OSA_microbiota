# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-09-06

# Last update: 2021-10-15

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
  
  
  # Rename levels for the factor variables OSAcat and t90cat (categories of sleep apnea)
  pheno[OSAcat=="no OSA", OSAcat:= "No Sleep Ap."]
  pheno[,OSAcat:=factor(OSAcat, levels = c("No Sleep Ap.", "Mild", "Moderate", "Severe"))]
  
  lev.t90cat <- levels(pheno[,t90cat])
  pheno[,t90cat:=factor(t90cat, levels=lev.t90cat, labels = c("0","T1","T2","T3"))]

  # Import results 
  res.m2 <- fread(paste0(input1,"cor2_all.var_mgs.tsv"))
  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]

  # Select the relevant MGSs 
  mgs.bmi <- res.m2[exposure =="BMI" & q.value<.05,MGS] # MGS correlated to BMI
  res.ahi <- res.m2[exposure =="ahi" & q.value<.05,] %>% filter(!MGS %in% mgs.bmi)
  mgs.ahi <- res.ahi$MGS
  res.t90 <- res.m2[exposure =="t90" & q.value<.05,] %>% filter(!MGS %in% mgs.bmi)
  mgs.t90 <- res.t90$MGS
  
  
  
  clean.y.axis <- function(y){
  y <- gsub("Absicoccus_", "A. ", y)
  y <- gsub("Bifidobacterium_longum_subsp._longum", "B. longum subsp.\nlongum", y)
  y <- gsub("Blautia_","B. ",y)
  y <- gsub("Pediococcus_", "P. ", y)
  y <- gsub("AM42_11","", y)
  y <- gsub("Staphylococcus_","S. ", y)
  y <- gsub("_sp"," sp", y)
  }
  
  # Result list 
  res.list <- list(ahi = res.ahi, t90 = res.t90)

  res.list <- lapply(res.list, function(x){x[,rho:=round(cor.coefficient,3)]})
  
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
    dades <- copy(data)
    
    setDF(dades)
    
    dades <- dades[!is.na(dades[,group]),]
    
    n.gr <- table(dades[,group])
    
    dades <- dades %>% group_by(get(group)) %>% summarise_at(pa,function(x){round((sum(x)/length(x))*100,1)})
    names(dades)[1] <- group
    names(dades) <- gsub("pa_","",names(dades))
  
    r <- res.list[[exposure]][MGS==y.axis,rho]
    
    ggplot(dades) + geom_bar(aes_string(x=group, y=y.axis), stat = 'identity') + 
      ggtitle(gsub("____","\n",clean.y.axis(y.axis)), subtitle = paste0("\u03C1=",r)) + 
      scale_x_discrete(labels=names(n.gr)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = .5, face = 'bold',size=7),
          plot.subtitle = element_text(hjust = .5, face = 'bold', size=8),
          axis.text.x = element_text(size = 7, angle = 45,hjust = 1),
          axis.text.y = element_text(size = 5), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank()) 

  }

  # Creating the plots 
  if(dir.exists("bar_ahi")==F){dir.create("bar_ahi")}
  
  
  message("Plot OSAcat")
  plots.ahi <- lapply(mgs.ahi, 
                    bar.plot.prevalence.fun, 
                    exposure = "ahi",
                    group = "OSAcat",
                    data = pheno[valid.ahi=='yes',])

  names(plots.ahi) <- mgs.ahi
  saveRDS(plots.ahi, file = paste0(output.plot,'/bar_ahi/barplot_ahi.rds'))

 # for(i in 1:length(plots.ahi)){
  #    ggsave(file = paste0("bar_ahi_",names(plots.ahi[i]),".png"), plot=plots.ahi[[i]],
    #     path = './bar_ahi/')
   #   }

  if(dir.exists("bar_t90")==F){dir.create("bar_t90")}
  
  plots.t90 <- lapply(mgs.t90, 
                    bar.plot.prevalence.fun, 
                    group= "t90cat",
                    exposure = "t90",
                    data = pheno[valid.t90=='yes',])

  names(plots.t90) <- mgs.t90
  saveRDS(plots.t90, file = paste0(output.plot,'/bar_t90/barplot_t90.rds'))

  #for(i in 1:length(plots.t90)){
  
   # ggsave(file = paste0("bar_t90_",names(plots.t90[i]),".png"), plot=plots.t90[[i]],
     #    path = './bar_t90/')
#  }


 # source('Script8_minibarplots.R')
  
  
  
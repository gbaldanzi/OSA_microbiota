# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-09-06

# Last update: 2021-09-06

# To investigate non-linear associations between MGS and AHI


rm(list = ls())
pacman::p_load(data.table,tidyverse,Hmisc)

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

pheno[valid.ahi=="yes", OSAcat := cut(ahi, breaks = c(0,5,15,30,Inf), include.lowest = T) ]
pheno[valid.t90=="yes", t90cat := cut(t90, breaks = c(0,10,15,Inf), include.lowest = T) ]

# Box plots 

# Function for plots 
box.plot.quantiles.fun <- function(y.axis,exposure,group,data){
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(Hmisc)
  dades <- data
  
  setDF(dades)
  
  #dades$quartile <-  cut(dades[[exposure]],
   #                      breaks = quantile(dades[[exposure]],
    #                                       probs = seq(0,1,by=.25),
     #                                      na.rm=T),
      #                   include.lowest = T)
  #dades <- dades %>% group_by(quartile) %>% 
    #mutate(mean.y=mean(get(y.axis),na.rm=T))
  
  r <- res.list[[exposure]][MGS==y.axis,rho]
  
  ggplot(dades) + geom_boxplot(aes_string(x=group, y=y.axis)) + 
    ggtitle(y.axis, subtitle = r) + 
    xlab(toupper(exposure)) + ylab(y.axis) +
   theme_classic() +
    theme(plot.title = element_text(hjust = .5, face = 'bold'),
          plot.subtitle = element_text(hjust = .5, face = 'bold'),
          axis.text.x = element_text(size = 8)) 

  }

plots.ahi <- lapply(mgs.ahi, 
                    box.plot.quantiles.fun, 
                    exposure = "ahi", 
                    group = "OSAcat",
                    data = pheno[valid.ahi=="yes",])

names(plots.ahi) <- mgs.ahi

dir.create("box_ahi")

for(i in 1:length(plots.ahi)){
  
  ggsave(file = paste0("box_ahi_",names(plots.ahi[i]),".png"), plot=plots.ahi[[i]],
         path = './box_ahi/')
}


plots.t90 <- lapply(mgs.t90, 
                    box.plot.quantiles.fun, 
                    exposure = "t90", 
                    group= "t90cat",
                    data = pheno[valid.t90=="yes",])

names(plots.t90) <- mgs.t90

dir.create("box_t90")

for(i in 1:length(plots.t90)){
  
  ggsave(file = paste0("box_t90_",names(plots.t90[i]),".png"), plot=plots.t90[[i]],
         path = './box_t90/')
}



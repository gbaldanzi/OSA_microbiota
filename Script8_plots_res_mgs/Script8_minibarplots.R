# Creating the plots 

library(ggplot2)
output.plot = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

plots.ahi <- readRDS(paste0(output.plot,'/bar_ahi/barplot_ahi.rds'))

for(i in 1:length(plots.ahi)){
  ggsave(file = paste0("minibar_ahi_",names(plots.ahi[i]),".png"), plot=plots.ahi[[i]],
         path = './bar_ahi/', dpi = 60)
}


#
  plots.t90 <- readRDS(paste0(output.plot,'/bar_t90/barplot_t90.rds'))

  for(i in 1:length(plots.t90)){
  
  ggsave(file = paste0("minibar_t90_",names(plots.t90[i]),".png"), plot=plots.t90[[i]],
         path = './bar_t90/', dpi=60)
  }

  
  #
  plots.ahi <- readRDS(paste0(output.plot,'/bar_ahi/barplot_quant_ahi.rds'))

  
  for(i in 1:length(plots.ahi)){
  
        ggsave(file = paste0("minibar_quant_ahi_",names(plots.ahi[i]),".png"), 
               plot=plots.ahi[[i]],
               path = './bar_ahi/', 
               dpi=60)
  }


  #
  plots.t90 <- readRDS( paste0( output.plot ,'/bar_t90/barplot_quant_t90.rds'))

  
  for(i in 1:length(plots.t90)){
  
        ggsave(file = paste0("minibar_quant_t90_",names(plots.t90[i]),".png"), plot=plots.t90[[i]],
              path = './bar_t90/', dpi=60)
  }
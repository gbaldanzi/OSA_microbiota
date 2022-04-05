# Creating the plots side by side 

library(ggplot2)
library(cowplot)
library(patchwork)
library(gridExtra)
library(data.table)
output.plot = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"


 n.gr <- table(pheno[valid.ahi=='yes',OSAcat])
  #n.gr <- c(1851,899,295,130)
  #names(n.gr) <- c("No Sleep Apnea", "Mild", "Moderate", "Severe")

  #subtitle <-paste0(paste0(names(n.gr),'=',n.gr,"; "),collapse = "")

  plots.ahi <- readRDS(paste0(output.plot,'/bar_ahi/barplot_ahi.rds'))

  plots.aligned1 = plot_grid(plotlist=plots.ahi,nrow=1)

  title = ggdraw() + draw_label("Signature species associated to AHI", fontface = 'bold',size=10) +  
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
          plot.background = element_rect(fill="white",color=NA))
  
  #subtitle = ggdraw() + draw_label(subtitle, fontface = 'bold', size = 7) +  
    #theme(plot.margin = unit(c(0,0,0,0),"cm"))
  
  plots.aligned1 = plot_grid(title, plots.aligned1,ncol=1,rel_heights = c(0.08,1))


  ggsave(file = "mergedbarplots_mgs_ahi.png", plot=plots.aligned1, path = './bar_ahi/',width = 8.5,
         height = 2.5)

  #ggsave(file = "mergedbarplots_mgs_ahi.png", plot=plots.aligned1, path = './bar_ahi/')

#
plots.t90 <- readRDS(paste0(output.plot,'/bar_t90/barplot_t90.rds'))

  n.gr <- table(pheno[valid.t90=='yes',t90cat])
 

  #subtitle <-paste0(paste0(names(n.gr),'=',n.gr,"; "),collapse = "")

  plots.aligned2 = plot_grid(plotlist=plots.t90,nrow=4,ncol = 6)

  title = ggdraw() + draw_label("Signature species associated to T90", fontface = 'bold', size =10) +  
    theme(plot.margin = unit(c(0.4,0.2,0.4,0.2),"cm"),
          plot.background = element_rect(fill="white",color=NA))

  #subtitle = ggdraw() + draw_label(subtitle, fontface = 'bold', size = 7) +  
    #theme(plot.margin = unit(c(0,0,0,0),"cm"))

  plots.aligned2 = plot_grid(title,plots.aligned2,ncol=1,rel_heights = c(0.08,4))

   ggsave(file = "mergedbarplots_mgs_t90.png", plot=plots.aligned2, path = './bar_t90/', width = 8.5, height=9)
  
   # ggsave(file = "mergedbarplots_mgs_t90.pdf", plot=plots.aligned2, path = './bar_t90/')


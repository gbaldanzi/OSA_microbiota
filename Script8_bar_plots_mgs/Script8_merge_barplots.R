# Creating the plots side by side 

library(ggplot2)
library(cowplot)
library(patchwork)
library(gridExtra)
output.plot = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

  n.gr <- table(pheno[valid.ahi=='yes',OSAcat])

  subtitle <-paste0(paste0(names(n.gr),'=',n.gr,"; "),collapse = "")

  plots.ahi <- readRDS(paste0(output.plot,'/bar_ahi/barplot_ahi.rds'))

  plots.aligned = plot_grid(plotlist=plots.ahi,nrow=1)

  title = ggdraw() + draw_label("Signature species correlated to AHI", fontface = 'bold',size=10) +  
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))
  
  subtitle = ggdraw() + draw_label(subtitle, fontface = 'bold', size = 7) +  
    theme(plot.margin = unit(c(0,0,0,0),"cm"))
  
  plots.aligned = plot_grid(title, subtitle,plots.aligned,ncol=1,rel_heights = c(0.08,0.05,1))

plots.aligned


  ggsave(file = "mergedbarplots_mgs_ahi.png", plot=plots.aligned, path = './bar_ahi/',width = 8.5,
         height = 2.5)


#
plots.t90 <- readRDS(paste0(output.plot,'/bar_t90/barplot_t90.rds'))

n.gr <- table(pheno[valid.t90=='yes',t90cat])

subtitle <-paste0(paste0(names(n.gr),'=',n.gr,"; "),collapse = "")

  plots.aligned = plot_grid(plotlist=plots.t90,nrow=4,ncol = 6)

  title = ggdraw() + draw_label("Signature species correlated to T90", fontface = 'bold', size =10) +  
    theme(plot.margin = unit(c(0.4,0.2,0.4,0.2),"cm"))

  subtitle = ggdraw() + draw_label(subtitle, fontface = 'bold', size = 7) +  
    theme(plot.margin = unit(c(0,0,0,0),"cm"))

  plots.aligned = plot_grid(title, subtitle,plots.aligned,ncol=1,rel_heights = c(0.08,0.06,4))

  plots.aligned

  ggsave(file = "mergedbarplots_mgs_t90.png", plot=plots.aligned, path = './bar_t90/', width = 8.5, height=9.8)


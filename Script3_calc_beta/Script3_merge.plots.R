# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-01-05

# Script to merge PCoA plots for AHI and T90 severity groups. 

  library(cowplot)
  library(data.table)
  library(ggplot2)

  source("Script3_calc_beta/Script3_calc_bc_ahi_plots.R")
  source("Script3_calc_beta/Script3_calc_bc_t90_plots.R")
  source("Script3_calc_beta/Script3_calc_bc_odi_plots.R")

  pcoa.plot.merged <- plot_grid(NULL,NULL,NULL,p1,p2,p3,NULL,NULL,NULL,labels = c("","","","A","B","C"), label_size = 12, nrow=3,
                                label_y = 0.97, ncol=3, rel_heights = c(.5,1,.4))

  ggsave("PCoA_sleepapnea.png", plot = pcoa.plot.merged, device = "png",
       path = "/home/baldanzi/Sleep_apnea/Results/")

  ggsave("PCoA_sleepapnea.pdf", plot = pcoa.plot.merged, device = "pdf", 
       path = "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/")
  
  
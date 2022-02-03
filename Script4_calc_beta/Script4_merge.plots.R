# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-01-05

# Script to merge PCoA plots for AHI, T90, and ODI severity groups. 

  library(cowplot)
  library(data.table)
  library(ggplot2)

  source(paste0(input.script, "Script4_calc_bc_ahi_plots.R"))
  source(paste0(input.script, "Script4_calc_bc_t90_plots.R"))
  source(paste0(input.script, "Script4_calc_bc_odi_plots.R"))

  pcoa.plot.merged <- plot_grid(NULL,NULL,NULL,p1,p2,p3,NULL,NULL,NULL,labels = c("","","","A","B","C"), label_size = 12, nrow=3,
                                label_y = 0.97, ncol=3, rel_heights = c(.5,1,.4))

  ggsave("PCoA_sleepapnea.pdf", plot = pcoa.plot.merged, device = "pdf", 
       path = output.res)
  
  
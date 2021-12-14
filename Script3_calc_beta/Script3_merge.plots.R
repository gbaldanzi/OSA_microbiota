# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# Last update: 2021-12-09

# Script to merge PCoA plots for AHI and T90 severity groups. 

  library(cowplot)

  source("Script3_calc_beta/Script3_calc_bc_ahi_plots.R")
  source("Script3_calc_beta/Script3_calc_bc_t90_plots.R")

  pcoa.plot.merged <- plot_grid(NULL,NULL,p1,p2,NULL,NULL,labels = c("","","A","B"), label_size = 12, nrow=3,
                                label_y = 0.97, ncol=2, rel_heights = c(.3,1,.4))

  ggsave("PCoA_sleepapnea.png", plot = pcoa.plot.merged, device = "png",
       path = "/home/baldanzi/Sleep_apnea/Results/")

  ggsave("PCoA_sleepapnea.pdf", plot = pcoa.plot.merged, device = "pdf", 
       path = "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/")
  
  
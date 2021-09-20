# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-12

  message("Plots LOESS - AHI and MGS") 

  dades = copy(valid.ahi)
  
  plot_loess = function (var.y,var.x="ahi",data=dades) {
    ggplot(data, aes_string(x = var.x, y = var.y)) + 
      geom_point() + 
      stat_smooth(geom = "smooth", method = 'loess', span=.7, se=T, formula = y~x) + 
      ggtitle(var.y) + 
      ylab(var.y) + xlab(var.x)
    }
  
  myplots <- lapply(mgs.ahi, plot_loess, data = dades, var.x="ahi")
  
  pdf(paste0(output.plot,"plot_loess_mgs_ahi.pdf"))
  ggarrange(plotlist = myplots, ncol = 2, nrow = 2)
  dev.off()

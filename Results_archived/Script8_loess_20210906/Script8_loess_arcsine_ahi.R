# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-25

  message("Plots LOESS - AHI and MGS") 

  plot_loess = function (var.y,var.x="ahi",data= plotdata.ahi) {
    ggplot(data, aes_string(x = var.x, y = var.y)) + 
      geom_point() + 
      stat_smooth(geom = "smooth", method = 'loess', span=.7, se=T, formula = y~x) + 
      ggtitle(var.y) + 
      ylab(var.y) + xlab(var.x)
    }
  
  myplots <- lapply(mgs.ahi, plot_loess, data = plotdata.ahi, var.x="ahi")
  
  pdf(paste0(output.plot,"plot_loess_arcsine_mgs_ahi.pdf"))
  ggarrange(plotlist = myplots, ncol = 2, nrow = 2)
  dev.off()

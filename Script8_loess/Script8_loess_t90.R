# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-23

message("Plots LOESS - T90 and MGS") 

dades = copy(valid.t90)

plot_loess = function (var.y,var.x="t90",data=dades) {
  ggplot(data, aes_string(x = var.x, y = var.y)) + 
    geom_point() + 
    stat_smooth(geom = "smooth", method = 'loess', span=.7, se=T, formula = y~x) + 
    ggtitle(var.y) + 
    ylab(var.y) + xlab(var.x)
}

myplots <- lapply(mgs.t90, plot_loess, data = dades, var.x="t90")

pdf(paste0(output.plot,"plot_loess_mgs_t90.pdf"))
ggarrange(plotlist = myplots, ncol = 2, nrow = 2)
dev.off()
# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-12

message("Plots - AHI and MGS") 

# Log+1 transformation 
dades = copy(valid.ahi)
dades[, (mgs.ahi) := lapply(.SD, function(x){log(x+1)}), .SDcols = mgs.ahi]

for(i in 1:length(mgs.ahi)){
  p <- ggplot(valid.ahi, aes_string("ahi",mgs.ahi[i])) + 
      geom_point() + 
      stat_smooth(geom = "smooth", method = 'loess', span=.2, se=T, formula = y~x)

  }

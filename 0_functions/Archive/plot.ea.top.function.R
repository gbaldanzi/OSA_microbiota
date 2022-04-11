# Top and bottom function 

plot.ea.top = function(data, 
                       ea.res,
                       pathway.name,
                       cor.var.name = "cor",
                       MGS.var.name = "MGS",
                       enrich.var.list,
                       output.plot
){
  require(fgsea)
  require(data.table)
  
  dd <- data
  enrich.var.list
  fgseaRES <- ea.res
  
  mgs.cor = dd[[cor.var.name]]
  names(mgs.cor) = dd[[MGS.var.name]]
  
  stopifnot(length(exp)==1)
  
 
  topPathways <- as.vector(fgseaRES[ES > 0][head(order(-pval),10), get(pathway.name)])
  g = plotGseaTable(enrich.var.list[topPathways], mgs.cor, fgseaRES, 
                gseaParam=0.5, render=F)
  return(g)
  }
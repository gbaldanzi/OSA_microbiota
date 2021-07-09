# Script 7 - Functions 

MGS.Enrich.Analysis <- function(data, 
                            cor.var.name = "cor",
                            MGS.var.name = "MGS",
                            enrich.var.list) {
  
    require(fgsea)
    require(data.table)
    dd <- data
    paths <- enrich.var.list
    if(any(class(dd) %in% "data.table")){setDF(dd)}
  
       # Naming the pvalue variable for the fgsea analysis 
    mgs.cor = dd[[cor.var.name]]
    names(mgs.cor) = dd[[MGS.var.name]]
  
    #temp=strsplit(names(mgs.pval), "____")
    #temp=matrix(unlist(temp), ncol=2, byrow=TRUE)
    #names(mgs.pval)=temp[,2]
    
    set.seed(123)
    res<- fgsea(pathways = paths, stats = rank(mgs.cor), eps = 0, scoreType = "pos",nPermSimple = 100000)
    setnames(res,"padj","q.value")
    res[,exposure:=unique(dd$exposure)]
    rm(mgs.cor)
    rm(dd)
    rm(paths)
    return(res)
}

  
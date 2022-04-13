# Functions for enrichment analysis 

MGS.Enrich.Analysis <- function(data, 
                            p.value.name = "p.value",
                            MGS.var.name = "mgs",
                            enrich.var.list,
                            direction = "positive") {
  
    dd <- data
    paths <- enrich.var.list
    if(any(class(dd) %in% "data.table")){setDF(dd)}
  
    # Naming the pvalue variable for the fgsea analysis 
    mgs.pval = dd[[p.value.name]]
    names(mgs.pval) = dd[[MGS.var.name]]
  
    # Stratified analysis by the direction of the correlation coefficients
    if(direction == "positive"){
      mgs.pval = mgs.pval[which(names(mgs.pval) %in% dd[dd$rho>0,MGS.var.name])]  
    }
    
    if(direction == "negative"){
      mgs.pval = mgs.pval[which(names(mgs.pval) %in% dd[dd$rho<0,MGS.var.name])]  
    }
    
    set.seed(123)
    res<- fgsea(pathways = paths, stats = rank(-mgs.pval), eps = 0, scoreType = "pos")
    
    setnames(res,"padj","q.value")
    
    res[,exposure := unique(dd$exposure)]
    res[,direction := direction]
    
    return(res)
}


  cutlast <- function(char,n){
  l <- nchar(char)
  a <- l-n+1
  return(substr(char,a,l))
  }

  
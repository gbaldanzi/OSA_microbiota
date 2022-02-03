# Functions for enrichment analysis 

# Gabriel Baldanzi 

MGS.Enrich.Analysis <- function(data, 
                            p.value.name = "p.value",
                            cor.var.name = "cor",
                            MGS.var.name = "MGS",
                            enrich.var.list,
                            direction = "both",
                            maxSize=650) {
  
    require(fgsea)
    require(data.table)
    dd <- data
    paths <- enrich.var.list
    if(any(class(dd) %in% "data.table")){setDF(dd)}
  
    if(isFALSE(direction %in% c("both", "positive", "negative"))){stop("Direction must be character 'both','positive' or 'negative'")}
  
    # Naming the pvalue variable for the fgsea analysis 
    mgs.pval = dd[[p.value.name]]
    names(mgs.pval) = dd[[MGS.var.name]]
  
    if(direction == "positive"){
      mgs.pval = mgs.pval[which(names(mgs.pval) %in% dd[dd[[cor.var.name]]>0,MGS.var.name])]  
    }
    
    if(direction == "negative"){
      mgs.pval = mgs.pval[which(names(mgs.pval) %in% dd[dd[[cor.var.name]]<0,MGS.var.name])]  
    }
    
    #temp=strsplit(names(mgs.pval), "____")
    #temp=matrix(unlist(temp), ncol=2, byrow=TRUE)
    #names(mgs.pval)=temp[,2]
    
    set.seed(123)
    res<- fgsea(pathways = paths, stats = rank(-mgs.pval), eps = 0, scoreType = "pos",maxSize=maxSize,nPermSimple = 100000)
    setnames(res,"padj","q.value")
    res[,exposure:=unique(dd$exposure)]
    rm(mgs.pval)
    rm(dd)
    rm(paths)
    return(res)
}


  cutlast <- function(char,n){
  l <- nchar(char)
  a <- l-n+1
  return(substr(char,a,l))
  }

  
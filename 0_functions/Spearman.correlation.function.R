# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Function: Spearman correlation 

spearman.function = function(x1, x2, covari=NULL, data){
  require(ppcor)
  require(fastDummies)
  require(data.table)
  
  if(x2=="ahi"){
    dataset <- data[valid.ahi=="yes",]
  } else if(x2 %in% c("t90","odi")){
    dataset <- data[valid.t90=="yes"]
  } else {  dataset <- data }
  
  if(any(class(dataset) %in% "data.table")){setDF(dataset)}
  
  # Create dummy variables
  dataset2 <- as.data.frame(model.matrix(~ ., dataset[,c(x2,covari)]))
  
  if(length(x1)==1){dataset2[[x1]] <- dataset[match(row.names(dataset2),rownames(dataset)),x1]}
  
  if(length(x1)>1){
    dataset2 <- merge(dataset2, dataset[,x1], by=0, all.x=T)
  }
  

  # Final covariates names
  cov <- colnames(dataset2)[which(!colnames(dataset2) %in% c(x1, x2, "(Intercept)","Row.names"))]
  
  
  result=data.frame()
  
  for(i in 1:length(x1)){
    
    cc = complete.cases(dataset2[, c(x1[i], x2, cov)] )
    res=pcor.test(x=dataset2[cc, x1[i]],
                  y=dataset2[cc, x2],
                  z=dataset2[cc, cov],
                  method = "spearman")
    
    temp <- data.frame(x=x1[i],exposure=x2,
                                rho=res$estimate,
                                p.value=res$p.value,
                                N=res$n,
                                method = res$Method,
                                model = deparse(substitute(covari)),
                                covariates=paste(covari,collapse = "+"))
    
    result <- rbind(result,temp)
  }
  
  
  if(length(x1)>1){
    result$q.value = p.adjust(result[,4], method="BH")
    if(x2 %in% c("ahi","t90","odi")){  #No need for rounding with the metabolites data. 
    result$q.value[result$q.value>=0.001] <- round(result$q.value[result$q.value>=0.001],3)
    }
    result <-  result[order(result$q.value),]
  } 
  
  return(result)
}

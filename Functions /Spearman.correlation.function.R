# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Function: Spearman correlation 

spearman.function = function(x1, x2, covari=NULL, data){
  if(any(class(data) %in% "data.table")){setDF(data)}
  numeric.covari = names(data[covari])[sapply(data[covari],is.numeric)]
  factor.covari =  names(data[covari])[sapply(data[covari],is.factor)]
  factor.covari = c(factor.covari,
                    names(data[covari])[sapply(data[covari],is.character)])
  #Exclude rows with incomplete observation
  if(any(is.na(data[,x1]))){stop("NA in x1")}
  if(any(is.na(data[,x2]))){stop("NA in x2")}
  
  notexclude = which(apply(data[c(numeric.covari,factor.covari)], 1, function(x){all(!is.na(x))}))
  temp.data=data[notexclude,]
  
  if(length(factor.covari)>0){  #Dummies variables 
    temp.data = dummy_cols(temp.data, select_columns = factor.covari,
                           remove_most_frequent_dummy = T, 
                           remove_selected_columns = T)
    factor.covari = names(temp.data)[!names(temp.data) %in% names(data)]
    #covariates for final model 
    cov=c(numeric.covari,factor.covari)
  }
  #Partial Sperman correlation 
  result=data.frame(matrix(ncol=7, nrow=length(x1)))
  for(i in 1:length(x1)){
    res=pcor.test(temp.data[,x1[i]],temp.data[,x2],temp.data[,cov],
                  method = "spearman")
    a=data.frame(x=x1[i],y=x2,
                 estimate=res$estimate,
                 p.value=res$p.value,
                 N=res$n,
                 method = res$Method,
                 covariates=paste(covari,collapse = "+"))
    result[i,]=a
  }
  if(length(x1)>1){
    result$q.value = p.adjust(result[,4], method="BH")
  }
  return(result)
}
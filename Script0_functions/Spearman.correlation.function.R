# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Function: Spearman correlation 

spearman.function = function(x1, x2, covari=NULL, data){
  
  dataset <- data
  
  if(any(class(dataset) %in% "dataset.table")){setDF(dataset)}
  numeric.covari = names(dataset[covari])[sapply(dataset[covari],is.numeric)]
  factor.covari =  names(dataset[covari])[sapply(dataset[covari],is.factor)]
  factor.covari = c(factor.covari,
                    names(dataset[covari])[sapply(dataset[covari],is.character)])
  #Exclude rows with incomplete observation
  if(any(is.na(dataset[,x1]))){stop("NA in x1")}
  if(any(is.na(dataset[,x2]))){stop("NA in x2")}
  
  notexclude = which(apply(dataset[c(numeric.covari,factor.covari)], 1, function(x){all(!is.na(x))}))
  temp.dataset=dataset[notexclude,]
  
  if(length(factor.covari)>0){  #Dummies variables 
    temp.dataset = dummy_cols(temp.dataset, select_columns = factor.covari,
                           remove_most_frequent_dummy = T, 
                           remove_selected_columns = T)
    factor.covari = names(temp.dataset)[!names(temp.dataset) %in% names(dataset)]
    #covariates for final model 
    cov=c(numeric.covari,factor.covari)
  }
  #Partial Sperman correlation 
  result=dataset.frame(matrix(ncol=7, nrow=length(x1)))
  for(i in 1:length(x1)){
    res=pcor.test(temp.dataset[,x1[i]],temp.dataset[,x2],temp.dataset[,cov],
                  method = "spearman")
    a=dataset.frame(x=x1[i],y=x2,
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
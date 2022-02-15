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
  numeric.covari = names(dataset[covari])[sapply(dataset[covari],is.numeric)]
  factor.covari =  names(dataset[covari])[sapply(dataset[covari],is.factor)]
  factor.covari = c(factor.covari,
                    names(dataset[covari])[sapply(dataset[covari],is.character)])
  #Exclude rows with incomplete observation
  if(any(is.na(dataset[,x1]))){stop("NA in x1")}
  if(any(is.na(dataset[,x2]))){stop("NA in x2")}
  
  notexclude <-  complete.cases(dataset[,covari])
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
  result=data.frame()
  
  for(i in 1:length(x1)){
    
    res=pcor.test(temp.dataset[,x1[i]],temp.dataset[,x2],temp.dataset[,cov],
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
    result$q.value[result$q.value>=0.001] <- round(result$q.value[result$q.value>=0.001],3)
    result <-  result[order(result$q.value),]
  }
  
  return(result)
}
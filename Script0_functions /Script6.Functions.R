# Script for functions used at the Script6

# Function for cleaning the Spearman correlation results 
Clean.Correlation.Results <- function(res){
  require(data.table)
  if(!any(class(res) %in% "data.table")){setDT(res)}
     
     if(any(names(res) %in% "cor.coefficient")){
       setnames(res,"cor.coefficient","cor.")
     }
     
     a = c("MGS","exposure","cor.","p.value","q.value","N")   
     res = res[,a,with=F]
     
     res[,cor.:=round(cor.,3)]
     res[q.value>=0.001, q.value:=round(q.value, digits = 3)]
     res[p.value>=0.001, p.value:=round(p.value, digits = 3)]
     
     return(res)
     
}


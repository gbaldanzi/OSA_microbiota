# PERMANOVA function for beta-diversity 

# The PERMANOVA function 
PermanovaFunction = function(outcome, exposure, covari, data, distance = "bray", nodes = 1){
  dades <- data
  require(vegan)
  require(parallel)
  require(data.table)
  if(!any(class(dades) %in% "data.table")){setDT(dades)}
  if(nrow(dades)==0){stop("Data has nrow == 0")}
  if(ncol(dades)==0){stop("Data has ncol == 0")}
  if(length(exposure) != 1){stop("exposure should have length == 1 ")}
  
  # Permanova from vegan package cannot handle missing information 
  a = complete.cases(dades[,covari,with=F])
  
  if(length(covari) > 0){
    covariates = paste(covari,collapse = "+")
    formula = as.formula(paste0(outcome,'[a,] ~',exposure,'+',covariates))
  }
  
  if(length(covari) == 0){
    formula = as.formula(paste0(outcome,'~',exposure,))
  }
  
  result <- adonis2(formula,dades[a,], by="margin", method = distance ,permutations = 9999, parallel = nodes)
  variables = rownames(result)
  result=as.data.frame(result)
  result[1,"N"]=sum(a)
  result$variables = variables 
  result = result[,c("variables",names(result[,-ncol(result)]))]
  return(result)
}


# PERMANOVA parallel function 

Permanova.parallel.FUN <- function(outcome,exposure,model,data,nod=16){
  
  expo_p <<-  exposure
  
  stopifnot(class(outcome) == "character")
  stopifnot(class(exposure) == "character")
  
  stopifnot(any(class(get(outcome)) %in%  "matrix"))

  set.seed(123)
  nod=nod   # Number of workers to be used 
  cl = makeCluster(nod)
  clusterExport(cl, varlist = c(outcome,"expo_p",
                                deparse(substitute(data)),
                                deparse(substitute(model))))
  clusterEvalQ(cl, library(vegan))
  clusterEvalQ(cl, library(data.table))
  res = PermanovaFunction(outcome = outcome, exposure = exposure, 
                          covari = model, data = data, nodes = nod)
  stopCluster(cl)
  
  return(res)
  
}

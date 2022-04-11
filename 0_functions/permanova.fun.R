# PERMANOVA function for beta-diversity 

# The PERMANOVA function 
PermanovaFunction = function(y, exposure, covari, data, nodes = 1){
  dades <- data

  require(vegan)
  require(parallel)
  require(data.table)
  if(!any(class(dades) %in% "data.table")){setDT(dades)}
  if(length(exposure) != 1){stop("exposure should have length == 1 ")}
  
  # adonis2 from vegan package does not handle missing information 
  dades <-  dades[complete.cases(dades[,covari,with=F]), ]
  
  # Making sure that y and dades have the same observations 
  outcome.matrix <- y[dades$SCAPISid,dades$SCAPISid]
  
  covariates = paste(covari,collapse = "+")
  formula = as.formula(paste0('outcome.matrix ~',exposure,'+',covariates))
  
  result <- adonis2(formula,dades, by="margin", method = "bray" ,permutations = 9999, parallel = nodes)
  variables = rownames(result)
  result=as.data.frame(result)
  result[1,"N"] <- nrow(dades)
  result$variables = variables 
  result = result[,c("variables",names(result[,-ncol(result)]))]
  return(result)
}


# PERMANOVA parallel function 

Permanova.parallel.FUN <- function(y,exposure,model,data,nod=16){
  
  expo_p <<-  exposure
  
  stopifnot(class(exposure) == "character")
  
  stopifnot(any(class(y) %in%  "matrix"))

  set.seed(123)
  nod=nod   # Number of workers to be used 
  cl = makeCluster(nod)
  clusterExport(cl, varlist = c(deparse(substitute(y)),
                                "expo_p",
                                deparse(substitute(data)),
                                deparse(substitute(model))))
  clusterEvalQ(cl, library(vegan))
  clusterEvalQ(cl, library(data.table))
  res = PermanovaFunction(y = y, exposure = exposure, 
                          covari = model, data = data, nodes = nod)
  stopCluster(cl)
  
  return(res)
  
}

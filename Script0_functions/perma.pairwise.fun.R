# PERMANOVA function for beta-diversity 


# The PERMANOVA function 
PermanovaFunction = function(outcome, group1, group2, covari, data, distance = "bray", nodes = 1 , group_var= group_var){
  require(vegan)
  require(parallel)
  require(data.table)
  
  dades <- data
  
  if(!any(class(dades) %in% "data.table")){setDT(dades)}
  if(nrow(dades)==0){stop("Data has nrow == 0")}
  if(ncol(dades)==0){stop("Data has ncol == 0")}
  if(!any(class(get(outcome)) %in% c("matrix"))){stop("outcome should be a matrix ")}
  if(length(group1) != 1){stop("group1 should have length == 1 ")}
  if(length(group2) != 1){stop("group2 should have length == 1 ")}
  if(class(group1)!="character"){stop("group1 should be a character")}
  if(class(group2)!="character"){stop("group1 should be a character")}
  
  # Restricting BC matrix to participants we want to compare
  group.id <- dades[get(group_var) %in% c(group1,group2),SCAPISid]
  dades <- dades[get(group_var) %in% c(group1,group2),]
  
  bc <- BC[rownames(BC) %in% group.id, colnames(BC) %in% group.id]
  
  outcome <- "bc"
  
  # Creating a new exposure variable 
  exposure <- paste(group1,group2, sep = "_")
  dades[,paste(group1,group2, sep = "_"):=factor(get(group_var))]
  
  # Permanova from vegan package cannot handle missing information 
  
  if(length(covari) > 0){
    a = complete.cases(dades[,covari,with=F])
    covariates = paste(covari,collapse = "+")
    formula = as.formula(paste0(outcome,'[a,] ~',exposure,'+',covariates))
  }
  
  if(length(covari) == 0){
    a = 1:nrow(dades)
    formula = as.formula(paste0(outcome,'~',exposure))
  }
  
  message("Running PERMANOVA, be patient")
  result <- adonis2(formula,dades[a,], by="margin", method = distance ,permutations = 9999, parallel = nodes)
  variables = rownames(result)
  result=as.data.frame(result)
  result[1,"N"]=sum(a)
  result$variables = variables 
  result = result[,c("variables",names(result[,-ncol(result)]))]
  return(result)
}


# Preparing the PERMANOVA


perma.fun <- function(outcome, group1, group2, covari, data, nodes = 1, group_var= group_var){
set.seed(123)
data.char <- deparse(substitute(data))
cl = makeCluster(nodes)
clusterExport(cl, varlist = c("BC",data.char,"model1","model2","model3"))
clusterEvalQ(cl, library(vegan))
clusterEvalQ(cl, library(data.table))
message(paste("PERMANOVA OSA categories -",group1,group2))
message(" ")
res = PermanovaFunction(outcome = outcome, group1 = group1, group2 = group2, covari = model3, data = data, nodes = nodes, group_var= group_var)
return(res)
stopCluster(cl)
}



# Function to produce a final table with main results, with FDR adjustment 

clean.res <- function(list){
  
  require(data.table)
  
  res = lapply(list, function(x){
    x[1,c("variables","Pr(>F)")]
  })
  
  t <- do.call(rbind,res)
  names(t) <- c("Comparison", "p.value")
  t$q.value <- p.adjust(t$p.value , method = "BH")
  return(t)
}


# Wrapper function to run pairwise permanova comparisons 

pairwise.perma.fun <- function(outcome = "BC", group_var = "OSAcat", covari = model3, data = pheno, nodes = 16){

  if(!length(unique(data[[group_var]]))>2) stop("Grouping variables (group_var) has less than 3 levels")

    list.res = vector(mode = "list",length = 0)  

  
    l_grp <- levels(data[[group_var]])

    g1 = l_grp[1] ; g2 = l_grp[-1]

    for(i in 1:length(g2)){
        message(paste(g1,g2[i],sep = "_"))
        list.res[[paste(g1,g2[i],sep = "_")]] <-  perma.fun(outcome = outcome,group_var= group_var, group1 = g1, group2 = g2[i], covari = covari, data = data, nodes = nodes)
    }

    g1 = g2[1] ; g2 = g2[-1]

    for(i in 1:length(g2)){
        message(paste(g1,g2[i],sep = "_"))
        list.res[[paste(g1,g2[i],sep = "_")]] <-  perma.fun(outcome = outcome, group_var= group_var,group1 = g1, group2 = g2[i], covari = covari, data = data, nodes = nodes)
    }

    if(length(l_grp)>3) {  
        g1 = g2[1] ; g2 = g2[-1]
  
        for(i in 1:length(g2)){
            message(paste(g1,g2[i],sep = "_"))
            list.res[[paste(g1,g2[i],sep = "_")]] <-  perma.fun(outcome = outcome, group_var= group_var,group1 = g1, group2 = g2[i], covari = covari, data = data, nodes = nodes)
        }
    }
    return(list.res)
}
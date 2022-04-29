# Project Sleep apnea and gut microbiota 

# Function to run Spearman's correlation and estimate the 
# the standard error via bootstrap


cor.boot <- function (x, y, z, data, nrep = 1000, conf.level=0.95, ncor=16) {
  
  tab <- data
  setDF(tab)
  
  tab <- data.frame(model.matrix(~.,data = tab[,c(x,y,z)]))
  
  tab <- tab[complete.cases(tab), -grep("Intercept",names(tab)) ]

  fit <- ppcor::pcor.test(tab[,x], tab[,y], tab[,-which(names(tab) %in% c(x,y))], 
                          method = "spearman")
  
  # Standard error 
    cor.fun <- function(data, ind) {
      ppcor::pcor.test(data[ind, 1], data[ind, 2], data[ind, 3:ncol(data)], 
           method = "spearman")$estimate
    }
    
    set.seed(7)

    simul <- boot::boot(tab, cor.fun, R = nrep, parallel = "multicore", ncpus = ncor)

    
    se <- sd(simul$t)
    
    # 95% Confidence interval 
    tri <- sort(na.omit(simul$t))
    if (any(!is.finite(tri))) {tri <- tri[-which(!is.finite(tri))]}
    int <- (1-conf.level)/2
    int.inf <- floor(nrep*int)
    int.sup <- ceiling(nrep*(1-int))
    
    conf.int <- paste0("[",round(tri[int.inf],3),",",round(tri[int.sup],3),"]")
    
    result <- data.frame(exposure = x, outcome = y ,
                         rho = fit$estimate,
                         se = se, conf.int = conf.int , 
                         p.value = fit$p.value,
                         N = fit$n,
                         covariates = paste(z,collapse = "+"))
    
    return(result)

}

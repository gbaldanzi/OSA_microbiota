
      
      lm.mgs <- function(mgs,phenotype="ahi",covariates = model2, data = valid.ahi){
        formula <- as.formula(paste(mgs ,"~",phenotype,"+",paste(covariates,collapse = "+")))
        temp.res <- lm(formula = formula, data = data)
        a <-  model.frame(temp.res)
        N <- nrow(a[complete.cases(a),])
        res <- data.frame(phenotype = phenotype,
                          mgs=mgs,
                          as.list(summary(temp.res)$coefficients[phenotype,]),
                          N = N)
        names(res) <- c("phenotype","mgs","beta","SE","t","p.value", "N")
        return(res)
      }


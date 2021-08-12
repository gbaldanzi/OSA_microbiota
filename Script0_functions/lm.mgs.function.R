
      
      lm.mgs <- function(mgs,phenotype="ahi",covariates = model2){
        formula <- as.formula(paste(mgs ,"~",phenotype,"+",paste(covariates,collapse = "+")))
        temp.res <- lm(formula = formula, data = valid.ahi)
        res <- data.frame(phenotype = phenotype,
                          mgs=mgs,
                          as.list(summary(temp.res)$coefficients[phenotype,]))
        names(res) <- c("phenotype","mgs","beta","SE","t","p.value")
        return(res)
      }


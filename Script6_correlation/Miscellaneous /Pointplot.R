# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-31

# Inferential Statistics 

# This code will investigate the association between MGS and AHI with and 
# without adjustment for BMI 

# Loading packages 
  pacman::p_load(data.table, ppcor, fastDummies, vegan)

# Output folders 
  input1 = '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/'
  output = "/home/baldanzi/Sleep_apnea/Results/"

# Importing data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

# Transforming two-level factor variables into numeric variables 
  dades = copy(pheno[valid.ahi=="yes",])
  a= c("Sex", "metformin","hypermed","dyslipmed","ppi")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

  # Removing rare MGS   
  # Calculating MGS prevalence 
  noms=grep("____",names(pheno),value=T)
  # presence-absence transformation: If a species is present, it becomes 1. If absent, becomes zero
  data_pa <- decostand(x = pheno[,noms,with=F], "pa")
  # calculate sum per species
  data_sum <- data.frame(prevalence=apply(data_pa, 2, sum),
                         percentage = apply(data_pa,2,sum)/nrow(pheno))
  data_sum$MGS = rownames(data_sum)
  a = data_sum$MGS[data_sum$prevalence<5] #19 MGS are only present in less than 5 individuals
  #a = data_sum$MGS[data_sum$percentage<.01] #383 MGS are present in less than 1% of individuals
  
  dades <- dades[ , -a, with=F] 
  
    # Correlations 
  source(paste0(input1,"Spearman.correlation.function.R"))
  
  
# Correlation between AHI and MGS - Step1 ####

#Preparing exposure, outcomes, and covariates
exposure="ahi"
outcomes=grep("___",names(dades),value=T)

#Covariates 
# model 1 : adjust for age + sex + alcohol + smoking + plate + received 
model1 <-   c("age", "Sex", "plate","shannon")
model2 <-   c(model1, "BMI")

# Running correlation 
  
message("Correlation between MGS and AHI")
  res<-  lapply(list(model1,model2),
              spearman.function,x1=outcomes,x2=exposure,data = dades)

  res[[1]]$model <- "model1"
  res[[2]]$model <- "model2"

  res1 <- do.call(rbind,res)

  names(res1) = c("MGS", "exposure", "cor.coefficient", "p.value", 
               "N", "method", "covariates","q.value","model")

  #Scatter plot
  
  data.plot <- res1 %>% select(MGS,cor.coefficient,model) %>% 
    pivot_wider(id_cols = MGS, names_from = model, 
                values_from = cor.coefficient)
  
  data.plot$mgs.ahionly <- 0
  data.plot$mgs.ahionly[data.plot$MGS %in% mgs.ahionly] <- 1
  
  data.plot <-  data.plot[data.plot$MGS %in% mgs.m1,]
  
  
  p <- ggplot(data.plot, aes(x=model1, y = model2, color=mgs.ahionly)) + 
    geom_point() + geom_abline(intercept = 0, slope = 1)
  
  p
  ggsave('cor.coef_with_withoutBMI.jpeg')
  
  
  # Point plot 
  
  data.plot <- res1 %>% select(MGS,cor.coefficient,model,q.value) %>% 
    pivot_wider(id_cols = MGS, names_from = model, 
                values_from = c(cor.coefficient,q.value))
  
  data.plot$mgs.ahionly <- 0
  data.plot$mgs.ahionly[data.plot$MGS %in% mgs.ahionly] <- 1
  
  data.plot <-  data.plot[data.plot$MGS %in% mgs.m1,]
  
  
  p <- ggplot(data.plot) +
    geom_point(aes(x=order(-q.value_model1), y = cor.coefficient_model1, shape=as.factor(mgs.ahionly)),
               color="green") + 
    geom_point(aes(x=order(-q.value_model1), y = cor.coefficient_model2, shape=as.factor(mgs.ahionly)),
               color="orange")
    
  p
  ggsave('cor.coef_with_withoutBMI_pointplot.jpeg')

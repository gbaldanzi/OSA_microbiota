# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-26

# To investigate non-linear associations between MGS and AHI


rm(list = ls())
pacman::p_load(data.table,Hmisc, rms,tidyverse,quantreg)

# input and output folders 
input1 = "/home/baldanzi/Sleep_apnea/Results/"
input2 = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

#output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"
output.plot = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

# Import full data 
valid.ahi <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/validsleep_MGS.shannon_Upp.rds")
valid.t90 <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/valid.t90_MGS.shannon_Upp.rds")

# Import results 
res.ahi <- fread(paste0(input1,"cor2_ahi_mgs.tsv"))
res.t90 <- fread(paste0(input1,"cor2_t90_mgs.tsv"))
res.bmi <- fread(paste0(input1,"cor2_BMI_mgs.tsv"))

res.ahi[,q.value:=round(q.value,3)]
res.t90[,q.value:=round(q.value,3)]
res.bmi[,q.value:=round(q.value,3)]

# Select the relevant MGSs 
mgs.bmi <- res.bmi[q.value<.05,]
mgs.ahi <- res.ahi[q.value<.05 & !MGS %in% mgs.bmi$MGS,MGS]
mgs.t90 <- res.t90[q.value<.05 & !MGS %in% mgs.bmi$MGS,MGS]

# Result list 
res.list <- list(ahi = res.ahi[MGS %in% mgs.ahi,],
                 t90 = res.t90[MGS %in% mgs.t90,])

res.list <- lapply(res.list, function(x){x[,rho:=paste0("rho = ",round(cor.coefficient,3))]})

# Quantile regression 

# Function for plots 
quantile.splines.plot.fun <- function(y.axis,x.axis,data){
  
  dades <- data
  
  if(any(class(dades) %in% "data.table")){
    setDF(dades)
  }
  
  form <-  as.formula(paste0(y.axis,"~ rcs(",x.axis,",4)"))
  fitq75 <- rq(formula = form, tau = .75,  data = dades)
  fitq90 <- rq(formula = form, tau = .90,  data = dades)
  fitq95 <- rq(formula = form, tau = .95,  data = dades)
  fitq99 <- rq(formula = form, tau = .99,  data = dades)
  
  newdf <- NULL
  newdf <- as.data.frame(seq(min(dades[x.axis]), max(dades[x.axis]), by=1))
  names(newdf) <- x.axis
  
  #newdf$q75rq <- predict.rq(fitq75, newdata = newdf)
  #newdf$q90rq <- predict.rq(fitq90, newdata = newdf)
  #newdf$q95rq <- predict.rq(fitq95, newdata = newdf)

  dades$q75 <- predict.rq(fitq75)
  dades$q90 <- predict.rq(fitq90)
  dades$q95 <- predict.rq(fitq95)
  dades$q99 <- predict.rq(fitq99)

  dades <- dades %>% select(all_of(c(x.axis,y.axis)),q75,q90,q95,q99) %>% 
    pivot_longer(cols = q75:q99, names_to = "quantile", values_to = "value")
  
  r <- res.list[[x.axis]][MGS==y.axis,rho]

  ggplot(dades) + geom_point(aes_string(x=x.axis,y=y.axis), size=.8) + 
    geom_line(data = dades, aes(x=get(x.axis),y=value,group=quantile,color=quantile),size=1.3,alpha=.7) + 
    ggtitle(y.axis, subtitle = r) + xlab(toupper(x.axis)) + ylab(y.axis) +
   theme_classic() +
    theme(plot.title = element_text(hjust = .5, face = 'bold'),
          plot.subtitle = element_text(hjust = .5, face = 'bold')) 

#ggplot(dades) + geom_point(aes_string(x=x.axis,y=y.axis)) + 
  #geom_line(data = newdf, aes(x=get(x.axis),y=value,group=quantile,color=quantile)) + 
  #ggtitle(y.axis) + xlab(x.axis) + ylab(y.axis)

  }

plots.ahi <- lapply(mgs.ahi, 
                    quantile.splines.plot.fun, 
                    x.axis = "ahi", 
                    data = valid.ahi)

pdf(file = paste0(output.plot,"plots_quantile_ahi.pdf"))
plots.ahi
dev.off()

plots.t90 <- lapply(mgs.t90, 
                    quantile.splines.plot.fun, 
                    x.axis = "t90", 
                    data = valid.t90)

pdf(file = paste0(output.plot,"plots_quantile_t90.pdf"))
plots.t90
dev.off()


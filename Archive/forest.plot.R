# Random forest plot of estimates and CI 

library(data.table)
library(tidyverse)

  # Folder 
  results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'

  # Results
  res <- fread(paste0(results.folder,"cor_hb_stratified.tsv"))
  res <- cor_hb_stratified
  setDT(res)

  
  # Plot
  new.var.cols <- c("rho", "se", "conf.int", "p.value" ,"covariates")
  res.long <- melt(res, measure = lapply(new.var.cols, function(x) paste0(x,c("_low","_high"))), 
                   value.name = new.var.cols, variable.name = "Hbgroup")
  res.long[, Hbgroup := factor(Hbgroup, levels = c(1,2), labels = c("low","high"))]
  
  # Extract first number from the confidence interval 
  res.long[, LCI := as.numeric(str_match(conf.int, "^\\[(-?[\\d\\.]+)")[,2])]
  
  # Extract second number from the confidence interval
  res.long[, UCI := as.numeric(str_extract(conf.int, "(?<=,).*[0-9]"))]

  res.long %>% filter(exposure == "t90") %>% ggplot(aes(x=rho, y= Hbgroup, xmin=LCI, xmax=UCI, col=Hbgroup)) +
    geom_pointrange() + geom_vline(xintercept = 0, lty="dashed", size=.5) +
    geom_point(size = 1) + 
    geom_linerange(aes(xmin=LCI, xmax=UCI)) +
    facet_wrap(~ outcome, nrow = length(unique(res.long$outcome)), strip.position = "left") + ylab("Species") +
    xlab("Spearman's coefficient") + 
    theme(strip.text.y.left = element_text(hjust = 0, vjust = 0.5, angle=0, face = "bold", size = 6),
          axis.text.y = element_blank(), axis.ticks.y = element_blank() )
  
  
  
  
  
              

  
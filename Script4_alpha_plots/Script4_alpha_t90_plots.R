# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# Lastest update: 2021-09-23

# This script will create plots and table describing alpha in participants with valid T90 measurement

# Loading packages
library(ggplot2)


# Import data
  pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")
  valid.t90 <- pheno[valid.t90=='yes',]
  
# Scatter plot: Shannon index against T90
  a = with(valid.t90, cor(t90, shannon, method = "spearman"))
  a = round(a,3)
  p1 = ggplot(data=valid.t90,aes(x=t90, y= shannon)) + 
    geom_point(color="lightskyblue") + ggtitle(paste0("Shannon index and T90 \nSpear. cor",a)) + 
    xlab("T90") +
    ylab("Shannon index") +
    geom_smooth(method='lm', alpha=.8) +
    theme_light() + 
    theme(plot.title = element_text(hjust=.5, size=14, face = "bold"),
          axis.title = element_text(size=14))
  
  ggsave("scatter.shannon.T90.png", plot = p1, device = "png", 
         path = "/home/baldanzi/Sleep_apnea/Descriptive/")


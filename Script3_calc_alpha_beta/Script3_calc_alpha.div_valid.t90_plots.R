# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics 

# This script will create plots and table describing alpha in participants with valid T90 measurement

# Loading packages


# Import data
  valid.t90 <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/valid.t90_MGS.shannon_Upp.rds")
  
# Scatter plot: Shannon index against T90
  a = with(valid.t90, cor(t90, shannon, method = "spearman"))
  a = round(a,3)
  p1 = ggplot(data=valid.t90,aes(x=t90, y= shannon)) + 
    geom_point(color="lightskyblue") + ggtitle(paste0("Shannon index and T90 \nSpear. cor",a)) + 
    xlab("T90") +
    ylab("Shannon index") +
    geom_smooth(method='lm', alpha=.8) +
    theme_light() + 
    theme(plot.title = element_text(hjust=.5, size=14),
          axis.title = element_text(size=14))
  
  ggsave("scatter.shannon.T90.png", plot = p1, device = "png", 
         path = "/home/baldanzi/Sleep_apnea/Descriptive/")
  
#Box plot
  res=with(valid.t90, kruskal.test(shannon ~ t90cat))
  p2 = ggplot(data=valid.t90) + 
       geom_boxplot(aes(x=t90cat, y= shannon, fill=t90cat)) + 
       ggtitle(paste0("Shannon index by T90 groups\nKruskal-Wallis: p.value=",
                   formatC(res$p.value,format= "e", digits=2)))+
       xlab("T90") +
       ylab("Shannon index") +
       theme(title = element_text(hjust=.5, size=14),
             axis.title = element_text(size=12))
  
  ggsave("box.shannon.T90.png", plot = p2, device = "png", 
         path = "/home/baldanzi/Sleep_apnea/Descriptive/")
  
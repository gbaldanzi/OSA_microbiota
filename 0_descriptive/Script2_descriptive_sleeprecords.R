# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Latest update: 2022-02-01

# Script to create descriptive plots of the sleep monitoring data 

  # Folder containing the dataset
  input <- "/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording/"

  # Folder to output the plots
  output.plot="/home/baldanzi/Sleep_apnea/Descriptive/"

  # Import data
  sleep <-  readRDS(paste0(input,"sleep.rds"))


  # Histogram of flow monitoring duration 
  p1 = ggplot(data = sleep, aes(x=fldeutv_min)) + geom_histogram() + 
    geom_vline(xintercept = 240, color = "red")  +
    scale_x_continuous(breaks = c(180, 240, 300, 420, 540, 660), labels = c("3h", "4h", "5h", "7h", "9h", "11h")) +
    ggtitle("Duration of flow monitoring") +  xlab("") +
    theme(axis.text.x = element_text(color = c("black", "red", "black", "black", "black", "black")))

  ggsave("hist.flow", plot = p1, device = "png", 
         path = output.plot)
  
  # Histogram of saturation monitoring duration 
  p2 = ggplot(data = sleep, aes(x=spo2utv_min)) + geom_histogram() + 
    geom_vline(xintercept = 240, color = "red")  +
    scale_x_continuous(breaks = c(180, 240, 300, 420, 540, 660), labels = c("3h", "4h", "5h", "7h", "9h", "11h")) +
    ggtitle("Duration of saturation monitoring") +  xlab("") +
    theme(axis.text.x = element_text(color = c("black", "red", "black", "black", "black", "black")))

  ggsave("hist.saturation", plot = p2, device = "png", 
         path = output.plot)

  # Histogram AHI
  p3 = ggplot(data = sleep[valid.ahi=='yes',], aes(x=ahi)) + geom_histogram()  +
    ggtitle("Hist Apnea-hypoapnea index") +  xlab("") + 
    geom_vline(xintercept = 5, color = "black", linetype = "twodash") + 
    geom_vline(xintercept = 15, color = "black", linetype = "twodash") + 
    geom_vline(xintercept = 30, color = "black", linetype = "twodash")

  ggsave("hist.ahi.png",plot = p3, device = "png", 
         path = output.plot)
  
  # Histogram ODI
  p4 = ggplot(data = sleep[valid.t90=='yes',], aes(x=odi)) + 
    geom_histogram(color="black", fill="lightskyblue2")  +
    ggtitle("Hist. Oxygen desaturation index") +  xlab("") +
    theme_light() +
    theme(plot.title= element_text(hjust = 0.5, size = 16,face = "bold"))

  ggsave("hist.odi.png",plot = p4, device = "png", 
       path = output.plot)

  # Histogram T90 
  p5 = ggplot(data = sleep[valid.t90=='yes',], aes(x=t90)) + 
    geom_histogram(color="black", fill="lightskyblue2") +
    geom_density(alpha=.2, fill="#FF6666") +
    ggtitle("Hist - T90%") +  xlab("") +
    theme_light()+
    theme(plot.title= element_text(hjust = 0.5, size = 16,face = "bold"))

  ggsave("hist.sat90.png",plot = p5, device = "png", 
       path = output.plot)

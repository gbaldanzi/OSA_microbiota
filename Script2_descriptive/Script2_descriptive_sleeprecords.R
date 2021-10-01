# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Latest update: 2021-09-22

# Script to create descriptive plots of the sleep recording 

output.plot="/home/baldanzi/Sleep_apnea/Descriptive/"

# Import data
sleep <-  readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/sleep_recording/sleep.rds")


#histogram of flow monitoring duration 
p1 = ggplot(data = sleep, aes(x=fldeutv_min)) + geom_histogram() + geom_vline(xintercept = 240, color = "red")  +
  scale_x_continuous(breaks = c(180, 240, 300, 420, 540, 660), labels = c("3h", "4h", "5h", "7h", "9h", "11h")) +
  ggtitle("Duration of flow monitoring") +  xlab("") +
  theme(axis.text.x = element_text(color = c("black", "red", "black", "black", "black", "black")))

#histogram of saturation monitoring duration 
p2 = ggplot(data = sleep, aes(x=spo2utv_min)) + geom_histogram() + geom_vline(xintercept = 240, color = "red")  +
  scale_x_continuous(breaks = c(180, 240, 300, 420, 540, 660), labels = c("3h", "4h", "5h", "7h", "9h", "11h")) +
  ggtitle("Duration of saturation monitoring") +  xlab("") +
  theme(axis.text.x = element_text(color = c("black", "red", "black", "black", "black", "black")))


# Histogram AHI
p3 = ggplot(data = sleep[valid.ahi=='yes',], aes(x=ahi)) + geom_histogram()  +
  ggtitle("Hist Apnea-hypoapnea index") +  xlab("") + 
  geom_vline(xintercept = 5, color = "black", linetype = "twodash") + 
  geom_vline(xintercept = 15, color = "black", linetype = "twodash") + 
  geom_vline(xintercept = 30, color = "black", linetype = "twodash")

# Histogram ODI
p4 = ggplot(data = sleep[valid.t90=='yes',], aes(x=odi)) + geom_histogram()  +
  ggtitle("Hist Oxygen desaturation index") +  xlab("")

# Histogram Sat90% 
p5 = ggplot(data = sleep[valid.t90=='yes',], aes(x=sat90)) + 
  geom_histogram(color="black", fill="lightskyblue2") +
  geom_density(alpha=.2, fill="#FF6666") +
  ggtitle("Hist - T90%") +  xlab("") +
  theme_light()+
  theme(plot.title= element_text(hjust = 0.5, size = 16,face = "bold"))

ggsave("hist.flow", plot = p1, device = "png", 
       path = )

ggsave("hist.saturation", plot = p2, device = "png", 
       path = output.plot)

ggsave("hist.ahi.png",plot = p3, device = "png", 
       path = output.plot)

ggsave("hist.odi.png",plot = p4, device = "png", 
       path = output.plot)

ggsave("hist.sat90.png",plot = p5, device = "png", 
       path = output.plot)

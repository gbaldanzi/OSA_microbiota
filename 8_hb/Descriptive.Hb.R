# Relation of Hb/Ht with OSA and impact on correlations with GM.

# 2022-04-20

library(data.table)
library(tidyverse)
library(cowplot)

# Folders 
results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
work <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/work/'


# load functions
source("0_functions/Spearman.correlation.function.R")

# Import data
pheno <- readRDS(paste0(work,"pheno_sleep_mgs_shannon.rds"))


# Hb by group OSA
pheno %>% group_by(OSAcat) %>% summarise(mean = mean(HBFormattedResult,na.rm=T),
                                         median = median(HBFormattedResult, na.rm=T))

pheno %>% group_by(t90cat) %>% summarise(mean = mean(HBFormattedResult,na.rm=T),
                                         median = median(HBFormattedResult, na.rm=T))

pheno %>% group_by(odicat) %>% summarise(mean = mean(HBFormattedResult,na.rm=T),
                                         median = median(HBFormattedResult, na.rm=T))

fit1 <- summary(aov(formula = HBFormattedResult ~ OSAcat + age + Sex, data = pheno))

p1 <- ggplot(pheno, aes(x=OSAcat,y=HBFormattedResult)) + ylab("Hb") +
  geom_boxplot() + ggtitle("AHI categories",
                           subtitle = paste("ANOVA p-value =",
                                            formatC(fit1[[1]][["Pr(>F)"]][1], 2, format='e')))

fit2 <- summary(aov(formula = HBFormattedResult ~ t90cat + age + Sex, data = pheno))
                           
p2 <- ggplot(pheno, aes(x=t90cat,y=HBFormattedResult)) + ylab("Hb") +
      geom_boxplot() + ggtitle("T90 categories",
                                subtitle = paste("ANOVA p-value =",
                                                 formatC(fit2[[1]][["Pr(>F)"]][1], 2, format='e')))

fit3 <- summary(aov(formula = HBFormattedResult ~ odicat + age + Sex, data = pheno))

p3 <- ggplot(pheno, aes(x=odicat,y=HBFormattedResult)) + ylab("Hb") +
  geom_boxplot() + ggtitle("ODI categories",
                           subtitle = paste("ANOVA p-value =",
                                            formatC(fit3[[1]][["Pr(>F)"]][1], 2, format='e')))
merged.p <- cowplot::plot_grid(p1,p2,p3)

ggsave(filename="HB_by_group.png",merged.p)
 
                                                


pheno[HBFormattedResult < median(HBFormattedResult,na.rm=T),Hbgroups := 1]

pheno[HBFormattedResult >= median(HBFormattedResult,na.rm=T),Hbgroups := 2]


pheno[,.(AHI = median(ahi,na.rm=T), T90 = median(t90,na.rm=T)),by=Hbgroups]
      



main.model <-   c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
main.model.BMI <- c(main.model, "BMI")


outcomes  <-  readRDS(paste0(results.folder,'mgs.m1.rds'))
exposures <- c("t90")



res <- lapply(exposures,spearman.function, 
               x1=outcomes,
               covari = main.model.BMI,
               data = pheno)

res <- do.call(rbind,res)

setDT(res)
setnames(res,"x","MGS")


res1 <- lapply(exposures,spearman.function, 
                         x1=outcomes,
                         covari = main.model.BMI,
                         data = pheno[Hbgroups==1])

res1 <- do.call(rbind,res1)

setDT(res1)
setnames(res1,"x","MGS")
res1 <- res1[match(res$MGS,res1$MGS),]



res2 <- lapply(exposures,spearman.function, 
               x1=outcomes,
               covari = main.model.BMI,
               data = pheno[Hbgroups==2])

res2 <- do.call(rbind,res2)

setDT(res2)
setnames(res2,"x","MGS")
res2 <- res2[match(res$MGS,res2$MGS),]


table.res <- data.frame(MGS=res$MGS,
                        rho = round(res$rho,3),
                        p.value = formatC(res$p.value,2,format = "e"),
                        N = res$N,
                        rho_HB1 = round(res1$rho,3),
                        p.value_HB1 = formatC(res1$p.value,2,format = "e"),
                        N_HB1 = res1$N,
                        rho_HB2 = round(res2$rho,3),
                        p.value_HB2 = formatC(res2$p.value,2,format = "e"),
                        N_HB2 = res2$N)

table.res[order(-table.res$rho),]

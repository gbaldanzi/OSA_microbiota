# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# This code creates a Venn diagram of MGS associated with either BMI or T90


# Results saved at "/home/baldanzi/Sleep_apnea/Results/discordant_t90.tsv" and 
# "/home/baldanzi/Sleep_apnea/Results/discordant_bmi_t90.tsv"

  rm(list = ls())
# Loading packages 
  pacman::p_load(data.table,ggplot2,ggVennDiagram, tidyr)

# Input folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  
# Output plot folder 
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Import data
  setwd(input)
  cor.t90 <- fread(paste0(input,"discordant_t90.tsv"))
  cor.bmi <- fread(paste0(input,"discordant_bmi_t90.tsv"))
  
# Prepare data for diagram 
  dades.sig.m1 <-  list(
    t90 = cor.t90[model=="model1" & q.value<0.05, MGS],
    BMI = cor.bmi[model=="model1" & q.value<0.05, MGS]
    )
  
  dades.sig.m2 <-  list(
    t90 = cor.t90[model=="model2" & q.value<0.05, MGS],
    BMI = cor.bmi[model=="model2" & q.value<0.05, MGS]
    )
  
  dades.sig.m3 <-  list(
    t90 = cor.t90[model=="model3" & q.value<0.05, MGS],
    BMI = cor.bmi[model=="model3" & q.value<0.05, MGS]
  )
  
  dades.sig.SA <-  list(
    t90 = cor.t90[model=="SA" & q.value<0.05, MGS],
    BMI = cor.bmi[model=="SA" & q.value<0.05, MGS]
    )
  
# Create VennDiagram    
  sig.m1 <- ggVennDiagram(dades.sig.m1) +
    ggtitle("MGS correlated with either BMI or t90 in model 1") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_bmit90_m1.png",plot = sig.m1,device = "png", path=output.plot)
  
  sig.m2 <- ggVennDiagram(dades.sig.m2) +
    ggtitle("MGS correlated with either BMI or t90 in model 2") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_bmit90_m2.png",plot = sig.m2,device = "png", path=output.plot)
  
  sig.m3 <- ggVennDiagram(dades.sig.m3) +
    ggtitle("MGS correlated with either BMI or t90 in model 3") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_bmit90_m3.png",plot = sig.m3,device = "png", path=output.plot)
  
  sig.SA <- ggVennDiagram(dades.sig.SA) +
    ggtitle("MGS correlated with either BMI or t90 /nmodel 2 excluding medication users") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  ggsave("Venn_bmit90_sa.png",plot = sig.SA,device = "png", path=output.plot)
  
# Discordant MGS associated with t90 but not with BMI  ####
  dis.m1 <- dades.sig.m1$t90[!dades.sig.m1$t90 %in% dades.sig.m1$BMI]
  dis.m2 <- dades.sig.m2$t90[!dades.sig.m2$t90 %in% dades.sig.m2$BMI]
  dis.m3 <- dades.sig.m3$t90[!dades.sig.m3$t90 %in% dades.sig.m3$BMI]
  dis.SA <- dades.sig.SA$t90[!dades.sig.SA$t90 %in% dades.sig.SA$BMI]
  
  stats = c("MGS","exposure", "cor.coefficient", "p.value","q.value", "N")
  a = cor.t90[model=="model1" & MGS %in% dis.m1, stats , with=F]
  b = cor.bmi[model=="model1" & MGS %in% dis.m1, stats , with=F]
  dis.m1 = rbind(a,b)
  dis.m1 <- dis.m1[order(dis.m1$MGS),]
  
  stats = c("MGS","exposure", "cor.coefficient", "p.value","q.value", "N")
  a = cor.t90[model=="model2" & MGS %in% dis.m2, stats , with=F]
  b = cor.bmi[model=="model2" & MGS %in% dis.m2, stats , with=F]
  dis.m2 = rbind(a,b)
  dis.m2 <- dis.m2[order(MGS),]
  
  stats = c("MGS","exposure", "cor.coefficient", "p.value","q.value", "N")
  a = cor.t90[model=="model3" & MGS %in% dis.m3, stats , with=F]
  b = cor.bmi[model=="model3" & MGS %in% dis.m3, stats , with=F]
  dis.m3 = rbind(a,b)
  dis.m3 <- dis.m3[order(MGS),]

  stats = c("MGS","exposure", "cor.coefficient", "p.value","q.value", "N")
  a = cor.t90[model=="SA" & MGS %in% dis.SA, stats , with=F]
  b = cor.bmi[model=="SA" & MGS %in% dis.SA, stats , with=F]
  dis.SA= rbind(a,b)
  dis.SA <- dis.SA[order(MGS),]
  
  discordant = list(model1=dis.m1,
                    model2=dis.m2,
                    model3=dis.m3,
                    SA=dis.SA)
  
  saveRDS(discordant, file = paste0(input,"discordant_t90_bmi.rds"))
  
  
  
# Discordant MGS associated with both t90 and BMI but opposite directions #### 
  both.m1 <- dades.sig.m1$t90[dades.sig.m1$t90 %in% dades.sig.m1$BMI]
  both.m2 <- dades.sig.m2$t90[dades.sig.m2$t90 %in% dades.sig.m2$BMI]
  both.m3 <- dades.sig.m3$t90[dades.sig.m3$t90 %in% dades.sig.m3$BMI]
  both.SA <- dades.sig.SA$t90[dades.sig.SA$t90 %in% dades.sig.SA$BMI]
  
  a = cor.t90[model=="model1" & MGS %in% both.m1, stats , with=F]
  b = cor.bmi[model=="model1" & MGS %in% both.m1, stats , with=F]
  both.m1 = rbind(a,b)
  both.m1 <- both.m1[order(MGS)]
  
  a = cor.t90[model=="model2" & MGS %in% both.m2, stats , with=F]
  b = cor.bmi[model=="model2" & MGS %in% both.m2, stats , with=F]
  both.m2 = rbind(a,b)
  both.m2 <- both.m2[order(MGS)]
  
  a = cor.t90[model=="model3" & MGS %in% both.m3, stats , with=F]
  b = cor.bmi[model=="model3" & MGS %in% both.m3, stats , with=F]
  both.m3 = rbind(a,b)
  both.m3 <- both.m3[order(MGS)]
  
  a = cor.t90[model=="SA" & MGS %in% both.SA, stats , with=F]
  b = cor.bmi[model=="SA" & MGS %in% both.SA, stats , with=F]
  both.SA= rbind(a,b)
  both.SA <- both.SA[order(MGS)]
  
  a = c("MGS", "exposure", "cor.coefficient")
  wide.m1 <- spread(both.m1[,a,with=F], exposure, cor.coefficient)
  wide.m2 <- spread(both.m2[,a,with=F], exposure, cor.coefficient)
  wide.m3 <- spread(both.m3[,a,with=F], exposure, cor.coefficient)
  wide.SA <- spread(both.SA[,a,with=F], exposure, cor.coefficient)
  
  # Check if there are MGS with opposite directions for the t90 or BMI 
  wide.m1[t90>0 & BMI<0, MGS]  # zero
  wide.m1[t90<0 & BMI>0, MGS]  # zero 
  wide.m2[t90>0 & BMI<0, MGS]  # zero
  wide.m2[t90<0 & BMI>0, MGS]  # zero
  wide.m3[t90>0 & BMI<0, MGS]  # zero
  wide.m3[t90<0 & BMI>0, MGS]  # zero
  wide.SA[t90>0 & BMI<0, MGS]  # zero
  wide.SA[t90<0 & BMI>0, MGS]  # zero
  
  
  
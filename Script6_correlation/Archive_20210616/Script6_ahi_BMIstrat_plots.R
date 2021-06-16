# Script6_ahi_BMIstrat_plots

rm(list = ls())


pacman::p_load(tidyverse,data.table,DT)

  input <- "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
  cor.ahi <- fread("/home/baldanzi/Sleep_apnea/Results/cor_ahi_mgs.tsv")

# Assign object with BMI groups and models
  group <-  c("<25","25-30",">=30")
  mo <- unique(cor.ahi[,model])

# Create plots and tables for every BMIcat

  for(i in 1:length(group)){
      cor.ahi <- fread(paste0(input,"cor_ahi_mgs_bmi",group[i],".tsv"), sep=",")
      setnames(cor.ahi,'cor.coeficient','cor')
      for(j in 1:length(mo)){
      a = sum(cor.ahi[model==mo[j],q.value]<.05,na.rm=T)
      n = nrow(cor.ahi[model==mo[j], ])

      p = ggplot(cor.ahi[model==mo[j],],aes(y=-log(q.value),x=phylum,color=phylum))+
      geom_jitter() + geom_hline(yintercept = -log(0.05), lty=2, color="darkred")+
      ggtitle(label=paste0("Spearman cor. AHI and MGS - ",mo[j]),
            subtitle = paste0("N = ",n," . Correlations FDR p-value<.05: ",a)) + 
              xlab("Phyla") + labs(tag = paste0("BMI",group[i])) +
      theme(plot.title = element_text(size=14, face="bold", hjust=.5),
            axis.text.x = element_text(angle=45,hjust = 1),
            legend.position = "none")
      
      ggsave(p, filename = paste0("res.plot_ahi_",i,"_",mo[j],".png"), device = "png",
             path = output.plot)

      t = cor.ahi[model==mo[j],]%>% select(MGS, exposure, cor, p.value,q.value,N) %>% arrange(q.value)
      t$cor = round(t$cor,digits = 3)
      t[q.value>=0.001, q.value:=round(q.value, digits = 3)]
      t[p.value>=0.001, p.value:=round(q.value, digits = 3)]

      saveRDS(t, file = paste0(output.plot,"res.table_ahi_",group[i],"_",mo[j],".rds"))
      
      }
  }

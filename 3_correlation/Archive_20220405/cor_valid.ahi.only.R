# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Created: 2022-02-01

# Sensitivity analysis, including only participants with valid AHI when assessing 
# the association between species and T90/ODI/BMI

# basic and full model will be used.  

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggvenn)

# Input and output folders 
  input <-  "/home/baldanzi/Datasets/sleep_SCAPIS/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Import data 
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  pheno <- pheno[valid.ahi='yes',] # Keep only participants with valid AHI
  
  
# Outcome variables 
  outcomes=grep("___",names(pheno),value=T)
  
#Covariates 
  basic.model <-  c("age", "Sex", "Alkohol","smokestatus","plate","shannon")
  full.model <- c(basic.model,"metformin","hypermed","dyslipmed","ppi","Fibrer",
                  "Energi_kcal" ,"leisurePA", "educat","placebirth","visit.month")
  
# Partial Spearman Correlation function 
  source("Script0_fuctions/Spearman.correlation.function.R")

# Basic model correlation ####
  
# AHI - basic model
  res.ahi <- fread(res, file = paste0(output,"cor_all.var_mgs.tsv"))
  res.ahi <- res.ahi[exposure=="ahi",]
  setnames(res.ahi, "MGS", "x")
  setDF(res.ahi)
  
# T90 - basic model
  res.t90 <-   spearman.function(x1=outcomes,
                                 x2="t90",covari = basic.model,
                                 data = pheno)
  
# ODI - basic model
  res.odi <-   spearman.function(x1=outcomes,
                                 x2=odi,covari = basic.model,
                                 data = pheno)
  
# BMI - basic model
  res.bmi <-   spearman.function(x1=outcomes,
                                 x2="BMI",covari = basic.model,
                                 data = pheno)
  
# Result list
  cols <- c("x","exposure","q.value")
    
  basic.model.list <- list(AHI=res.ahi[,cols], 
                           T90 = res.t90[,cols], 
                           ODI = res.odi[,cols], 
                           BMI = res.bmi[,cols])
  
# Species associated with the phenotypes in basic model 
  res.table <- do.call(rbind,basic.model.list)
  species.basic.model <- unique(res.table[res.table$q.value<.05,"x"]) 
  
# Full model correlations ####
  
  full.model.list <-   lapply(c("ahi","t90","odi","BMI"),
                                spearman.function,
                                x1=species.basic.model,
                                covari = full.model,
                                data = pheno)
  names(full.model.list) <- c("AHI","T90","ODI","BMI")
  
  # Venn diagram Basic model 
  
  mgs.fdr = lapply(basic.model.list,function(df){df[x$q.value<0.05,"x"]})
  
  venn.m1 <- ggvenn(mgs.fdr,
                    fill_color = c("orange","cornflowerblue","grey83","green4"),
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .6, 
                    set_name_size = 4,
                    text_size = 3) +
    ggtitle("Basic model") +
    theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))

  # # Venn diagram full model 
  
  mgs.fdr = lapply(full.model.list,function(df){df[x$q.value<0.05,"x"]})
  
  venn.m2 <- ggvenn(mgs.fdr,
                    fill_color = c("orange","cornflowerblue","grey83","green4"),
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .6, 
                    set_name_size = 4,
                    text_size = 3) +
    ggtitle("Full model") +
    theme(plot.title = element_text(size=12, face = "bold", hjust = 0.5))
  
  # Merge Venn diagrams for the both models 
  
  venn <- plot_grid(venn.m1, venn.m2, labels = c("A","B"), label_size = 12, nrow=1,
                    label_y = .7)
  
  ggsave("Venn_sleepapnea_valid.ahi.only.png",plot = venn,device = "png", path=output.plot)
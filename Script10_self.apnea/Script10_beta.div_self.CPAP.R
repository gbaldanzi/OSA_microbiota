# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics continued:

# This script will calculate beta diversity with self-reported sleep apnea 

rm(list=ls())

# Loading packages
library(ape)
library(data.table)
library(vegan)
library(ggplot2)

pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")


  pheno[,self.CPAP:='NO']
  pheno[cqhe061=='CPAP',self.CPAP:='YES']
  
  sev.ahi.indiv <- pheno[OSAcat=="Severe",]
  mod.ahi.indiv <- pheno[OSAcat %in% c("Moderate","Severe"),]

# Bray-curtis dissimilarity #### 

  # Create dataset containing exclusively the MGS as columns/variables 
  a=grep("____",names(mod.ahi.indiv),value=T)
  MGSsev=sev.ahi.indiv[,a, with=F]
  MGSmod=mod.ahi.indiv[,a, with=F]

# Estimating the BC index - creates a matrix with the BC index between individual samples 
# MGS as relative abundances
  BCsev=as.matrix(vegdist(MGSsev,method="bray"))
  rownames(BCsev)=colnames(BCsev)=sev.ahi.indiv$SCAPISid
  fwrite(BCsev,'/home/baldanzi/Datasets/sleep_SCAPIS/sev.apnea.BCmatrix.csv',sep=",")
  
  BCmod=as.matrix(vegdist(MGSmod,method="bray"))
  rownames(BCmod)=colnames(BCmod)=mod.ahi.indiv$SCAPISid
  fwrite(BCmod,'/home/baldanzi/Datasets/sleep_SCAPIS/mod.apnea.BCmatrix.csv',sep=",")

  # Severe Apnea indiv ####
  
  # Principal coordinates analysis based on BC 
  pcoa.bray <- pcoa(BCsev)
  
  # Prepare data to PCoA plot 
  dat.plot=data.frame(pcoa.bray$vectors, CPAP = sev.ahi.indiv$self.CPAP)
  n = nrow(dat.plot)
  n1 = nrow(dat.plot[,dat.plot$CPAP=="NO"])
  n2 = nrow(dat.plot[,dat.plot$CPAP=='YES'])
  
  p=ggplot(dat.plot, aes(x=Axis.1,y=Axis.2, color=CPAP))+
    geom_point(size=1.1,aes(color=CPAP)) +
    stat_ellipse(type = "t", size=1.3) +
    scale_color_manual(values=c("gray84","blue3"),
                       labels=c(paste0("NO, n=",n1),
                                paste0("YES, n=",n2))) + 
    ggtitle(paste0("PCoA plot of Bray-curtis dissimilarity, N=",n),
            subtitle = "Severe measured apnea and self-reported CPAP treatment") +
    xlab(paste0("PCo1 \n (",round(100*pcoa.bray$values$Relative_eig[1],1),"% )")) +
    ylab(paste0("PCo2 \n (",round(100*pcoa.bray$values$Relative_eig[2],1),"% )")) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 0),
          plot.title = element_text(hjust = 0.5, face="bold", size=12),
          legend.text = element_text(size=10), 
          plot.subtitle = element_text(size=10, hjust = 0.5))
  
  ggsave("BC.sev.apnea.png", plot = p, device = "png",
        path = "/home/baldanzi/Sleep_apnea/Descriptive/")
  
  # Moderate and Severe Apnea Indiv ####
  
  # Principal coordinates analysis based on BC 
  pcoa.bray <- pcoa(BCmod)
  
  # Prepare data to PCoA plot 
  dat.plot=data.frame(pcoa.bray$vectors, CPAP = mod.ahi.indiv$self.CPAP)
  n = nrow(dat.plot)
  n1 = nrow(dat.plot[dat.plot$CPAP=='NO',])
  n2 = nrow(dat.plot[dat.plot$CPAP=='YES',])
  
  p2=ggplot(dat.plot, aes(x=Axis.1,y=Axis.2, color=CPAP))+
    geom_point(size=1.1,aes(color=CPAP)) +
    stat_ellipse(type = "t", size=1.3) +
    scale_color_manual(values=c("gray84","blue3"),
                       labels=c(paste0("NO, n=",n1),
                                paste0("YES, n=",n2))) +
    ggtitle(paste0("PCoA plot of Bray-curtis dissimilarity, N=",n),
            subtitle = "Moderate-severe measured apnea and self-reported CPAP treatment") +
    xlab(paste0("PCo1 \n (",round(100*pcoa.bray$values$Relative_eig[1],1),"% )")) +
    ylab(paste0("PCo2 \n (",round(100*pcoa.bray$values$Relative_eig[2],1),"% )")) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 0),
          plot.title = element_text(hjust = 0.5, face="bold", size=12),
          legend.text = element_text(size=10), 
          plot.subtitle = element_text(size=10, hjust = 0.5))
  
  ggsave("BC.mod.apnea.png", plot = p2, device = "png",
         path = "/home/baldanzi/Sleep_apnea/Descriptive/")
  
  
  
  # PERMANOVA ####
  # Loading packages 
  library(parallel)
  
  output = "/home/baldanzi/Sleep_apnea/Results/"
 
  source('permanova.fun.R')
  
  # Severe Apnea
  
  # Transforming two-level factor variables into numeric variables 
  dades = copy(sev.ahi.indiv)
  a= c("Sex")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]
  
  # Making sure that BC and dades have the same order of observations 
  dades = dades[match(rownames(BCsev),dades$SCAPISid),]
  
  # Outcome - character name (length=1) with matrix distance 
  outc = "BCsev"
  
  # Main Exposure - character name (length=1)
  expo = "self.CPAP"
  
  #Covariates 
  # model 1 : adjust for age + sex + alcohol + smoking + plate + received 
  model1 <-   c("age", "Sex")

  # PERMANOVA
  nod=4
  res = PermanovaFunction(outcome = outc, exposure = expo, 
                          covari = model1, data = dades, nodes = nod)
  fwrite(res, file = paste0(output,"permanova_sev.apnea_bc.tsv"), sep="\t")
  
  
  # Moderate and Severe Apnea
  
  # Transforming two-level factor variables into numeric variables 
  dades = copy(mod.ahi.indiv)
  a= c("Sex")
  dades[,(a):=as.data.frame(data.matrix(data.frame(unclass(dades[,a, with=F]))))]

  # Making sure that BC and dades have the same order of observations 
  dades = dades[match(rownames(BCmod),dades$SCAPISid),]
  
  # Outcome - character name (length=1) with matrix distance 
  outc = "BCmod"
  
  # Main Exposure - character name (length=1)
  expo = "self.CPAP"
  
  #Covariates 
  # model 1 : adjust for age + sex + alcohol + smoking + plate + received 
  model1 <-   c("age", "Sex")
  
  # PERMANOVA
  res = PermanovaFunction(outcome = outc, exposure = expo, 
                          covari = model1, data = dades, nodes = nod)
  fwrite(res, file = paste0(output,"permanova_mod.apnea_bc.tsv"), sep="\t")
  
  
  
  
  


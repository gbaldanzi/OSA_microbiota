# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi

# Last update: 2022-01-03

# This script will create a table showing the prevalence of each signature species 
# by category of AHI, T90 or ODI 

rm(list = ls())
pacman::p_load(data.table,tidyverse,Hmisc, vegan, cowplot)

# input and output folders 
input1 = "/home/baldanzi/Sleep_apnea/Results/"
input2 = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

#output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"
output.plot = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

# Import full data 
pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")


# Rename levels for the factor variables OSAcat and t90cat (categories of sleep apnea)
pheno[OSAcat=="no OSA", OSAcat:= "No Sleep Ap."]
pheno[,OSAcat:=factor(OSAcat, levels = c("No Sleep Ap.", "Mild", "Moderate", "Severe"))]

lev.t90cat <- levels(pheno[,t90cat])
pheno[,t90cat:=factor(t90cat, levels=lev.t90cat, labels = c("0","T1","T2","T3"))]

pheno[,odicat := as.factor( cut(pheno$odi,breaks = quantile(odi, probs = seq(0,1,by=.25), na.rm=T), 
                                include.lowest = T) )]


# Import results 
res.m2 <- fread(paste0(input1,"cor2_all.var_mgs.tsv"))
res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]
res.m2[,rho:=round(cor.coefficient,3)]

# Select the relevant MGSs 
mgs.fdr.m2 <- readRDS('/home/baldanzi/Sleep_apnea/Results/mgs.m2.rds')
mgs.fdr <- unique(do.call('c',mgs.fdr.m2))  

res.ahi <- res.m2[exposure =="ahi" & MGS %in% mgs.fdr.m2$mgs.fdr.ahi,] 
res.t90 <- res.m2[exposure =="t90" & MGS %in% mgs.fdr.m2$mgs.fdr.t90,] 
res.odi <- res.m2[exposure =="odi" & MGS %in% mgs.fdr.m2$mgs.fdr.odi,]



# Calculating prevalence 
pa <- paste0("pa_",mgs.fdr) # pa = present or absent 
pheno[,(pa):= lapply(.SD, decostand, method="pa"), .SDcols = mgs.fdr]

dades <- copy(pheno)
dades <- dades[valid.ahi =="yes",]

dades.ahi <- dades %>% group_by(OSAcat) %>% summarise_at(pa,function(x){round((sum(x)/length(x))*100,1)})
#dades.ahi[,1] <- group
dades.ahi[,1] <- c("a1","a2","a3","a4")
dades.ahi$exposure <- "AHI"
names(dades.ahi)[1] <- "Group"

dades <- copy(pheno)
dades <- dades[valid.t90 =="yes",]

dades.t90 <- dades %>% group_by(t90cat) %>% summarise_at(pa,function(x){round((sum(x)/length(x))*100,1)})
#dades.t90[,1] <- group
dades.t90[,1] <- c("t1","t2","t3","t4")
dades.t90$exposure <- "T90"
names(dades.t90)[1] <- "Group"

dades.odi <- dades %>% group_by(odicat) %>% summarise_at(pa,function(x){round((sum(x)/length(x))*100,1)})
#dades.odi[,1] <- group
dades.odi[,1] <- c("o1","o2","o3","o4")
dades.odi$exposure <- "ODI"
names(dades.odi)[1] <- "Group"

dades <- rbind(dades.ahi,dades.t90,dades.odi)

dades <- dades[,c("Group","exposure",pa)]

names(dades) <- gsub("pa_","",names(dades))

dades$exposure = factor(dades$exposure, levels = c("AHI","ODI", "T90"))


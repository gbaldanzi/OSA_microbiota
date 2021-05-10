# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

pacman::p_load(data.table,ppcor, fastDummies, vegan, ggplot2)

rm(list = ls())
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Script 3. Inferential statistics 

# Importing data
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
valid.ahi = fread("validsleep_MGS.shannon.BC_Upp.tsv",header=T, na.strings=c("", "NA")) 
setnames(valid.ahi, "pob", "placebirth")

# Correlation between AHI and Shannon ####

# Creating dummy variables 
valid.ahi.dummies = dummy_cols(valid.ahi, select_columns = c("plate","received","smokestatus"
                                                             ,"leisurePA","educat", "placebirth"), 
                               remove_most_frequent_dummy = T, remove_selected_columns = F)

# Transforming two-level factor variables into numeric variables 
a= c("Sex","ppi", "diabd","hypertension","dyslipidemia","diabmed","hypermed","dyslipmed")
valid.ahi.dummies[,(a):=as.data.frame(data.matrix(data.frame(unclass(valid.ahi.dummies[,a, with=F]))))]

# model 1 : adjust for age + sex + BMI + alcohol + smoking + plate + received 
model1 =  c("age", "Sex", "BMI", "Alkohol",
            grep("smokestatus_", names(valid.ahi.dummies), value=T),
            grep("plate_", names(valid.ahi.dummies), value=T), 
            grep("received_", names(valid.ahi.dummies), value=T))
# model 2 = model 1 + fiber intake + ppi + physical activity + education + country of birth 
model2 = c(model1, "Fibrer", "ppi", 
           grep("leisuraPA_", names(valid.ahi.dummies), value=T),
           grep("educat_", names(valid.ahi.dummies), value=T), 
           grep("placebirth_", names(valid.ahi.dummies), value=T))
# model 3 = model 2 + diabetes + hypertension + dyslipidemia, medication 
model3 = c(model2,"diabd","hypertension","dyslipidemia","diabmed","hypermed","dyslipmed")

models = list(model1=model1, model2=model2, model3=model3)

# Partial Spearman correlation loop for all 3 models ####
# pcor.test does not allow missing data
# Excluding individuals with missing information on covariates 
sp = data.frame(matrix(ncol=7, nrow=3))
for(i in 1:3){
notexclude = which(apply(valid.ahi.dummies[,models[[i]],with=F], 1, function(x){all(!is.na(x))}))
temp = pcor.test(valid.ahi.dummies[notexclude,shannon],
                valid.ahi.dummies[notexclude,ahi],
                valid.ahi.dummies[notexclude,models[[i]],with=F], 
                method = "spearman")
sp[i,] = c(as.data.frame(names(models[i])),temp)
}
names(sp) = c("model", names(temp))

# Saving results for correlation shannon and ahi 
fwrite(sp, file = paste0(output,"cor.shannon.ahi.tsv"), sep="\t")


# Bray-curtis dissimilarity #### 

# Investigating the association between BC and AHI 

# adjusting variables to factor format 
setwd("/home/baldanzi/Datasets/sleep_SCAPIS")
load("data_table1")
a = names(dat1[,-"SCAPISid"])
valid.ahi[,(a):=NULL]
valid.ahi = merge(dat1,valid.ahi, by = "SCAPISid", all=T, all.x=F, all.y=F)
valid.ahi$plate = as.factor(valid.ahi$plate)
valid.ahi$received = as.factor(valid.ahi$received)

# Importing BC matrix
BC = as.matrix(fread("OSA.BCmatrix.csv",header = T))  
rownames(BC) = colnames(BC)

# Permanova ####
set.seed(123)
perma1 <- adonis2(formula = BC ~ ahi + age + Sex + plate + received, 
                  data = valid.ahi, 
                  permutations = 1000, method = "bray", by = "margin")

perma2 <- adonis2(formula = BC ~ ahi + age + Sex + plate + received +
                    BMI + Alkohol + smokestatus, 
                  data = valid.ahi, 
                  permutations = 1000, method = "bray", by = "margin")

perma3 <- adonis2(formula = BC ~ ahi + age + Sex + plate + received +
                    BMI + Alkohol + smokestatus +
                    leisurePA + educat + placebirth + Fibrer, 
                  data = valid.ahi, 
                  permutations = 1000, method = "bray", by = "margin")

perma4 <- adonis2(formula = BC ~ ahi + age + Sex + plate + received +
                    BMI + Alkohol + smokestatus +
                    leisurePA + educat + placebirth + Fibrer +
                    diabd + hypertension + dyslipidemia + diabmed +
                    hypermed + dyslipmed , 
                  data = valid.ahi, 
                  permutations = 1000, method = "bray", by = "margin")

# Exporting PERMANOVA results 
perma1 = cbind(perma1, "model1")
perma2 = cbind(perma2, "model2")
perma3 = cbind(perma3, "model3")
perma4 = cbind(perma4, "model4")
fwrite(perma1, file = paste0(output,"permanova.tsv"), sep="\t")
fwrite(perma2, file = paste0(output,"permanova.tsv"), append =T,sep="\t")
fwrite(perma3, file = paste0(output,"permanova.tsv"), append =T,sep="\t")
fwrite(perma4, file = paste0(output,"permanova.tsv"), append =T,sep="\t")


# Correlation between AHI and MGS - Step1 ####
#only the model 1 will be used in this analysis 

sp2 = data.frame(matrix(ncol=7, nrow=1985))
j=0
for(i in grep("____H",names(valid.ahi.dummies))){
    j=j+1
  temp = pcor.test(valid.ahi.dummies[,(i), with=F],
                   valid.ahi.dummies[,ahi],
                   valid.ahi.dummies[,models[[1]],with=F], 
                   method = "spearman")
  sp2[j,] = cbind(names(valid.ahi.dummies[,(i),with=F]),temp)
}
names(sp2) = c("MGS", names(temp))

sp2 = sp2[order(sp2$p.value),]
fwrite(sp2, file = paste0(output,"MGScor_step1.tsv"),sep="\t")

# Filtering out rare species 
# Filter based on presence or absence 
# presence-absence transformation:
# If a species is present, it is transformed to 1
# If it is absent, it is transformed to zero
a = grep("____H",names(valid.ahi.dummies))
data_pa <- decostand(x = valid.ahi.dummies[,a,with=F], "pa")

# calculate sum per species
data_sum <- apply(data_pa, 2, sum)

p1 = ggplot(data = as.data.frame(data_sum), aes(x=data_sum)) + 
  geom_histogram(bins = 70)  +
  ggtitle("Histogram MGS prevalence") +  
  xlab(paste("Present in at least 100 indiv: ", 
            length(data_sum[data_sum>=100])," (",
            round(length(data_sum[data_sum>=100])/length(data_sum),2)*100,"%)\n",
            "Present in at least 50 indiv: ",
            length(data_sum[data_sum>=50])," (",
            round(length(data_sum[data_sum>=50])/length(data_sum),2)*100,"%)\n",
            sep="")) + 
  geom_vline(xintercept = 50, color = "black", linetype = "twodash") + 
  geom_vline(xintercept = 100, color = "black", linetype = "twodash")

ggsave("hist.MGS.prevalence",plot = p1, device = "png", 
       path = "/home/baldanzi/Sleep_apnea/Descriptive/")

#filter/remove rare species 
#remove species that do not contain at least 100 non-zero values
data_sum = as.data.frame(data_sum)
data_sum$MGS = rownames(data_sum)
a = data_sum$MGS[data_sum$data_sum<100]
v.ahi.dum.filt <- valid.ahi.dummies[ , -a, with=F]

# Correlation between AHI and MGS AFTER FILTERING - Step1 ####

#only the model 1 will be used in this analysis 

sp3 = data.frame(matrix(ncol=7, nrow=1194))
j=0

for(i in grep("____H",names(v.ahi.dum.filt))){
  j=j+1
  temp = pcor.test(v.ahi.dum.filt[,(i), with=F],
                   v.ahi.dum.filt[,ahi],
                   v.ahi.dum.filt[,models[[1]],with=F], 
                   method = "spearman")
  sp3[j,] = cbind(names(v.ahi.dum.filt[,(i),with=F]),temp)
}
names(sp3) = c("MGS", names(temp))

sp3 = sp3[order(sp3$p.value),]
fwrite(sp3, file = paste0(output,"MGScorfiltered_step1.tsv"),sep="\t")

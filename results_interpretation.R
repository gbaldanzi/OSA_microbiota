set.seed(101)

n <- 200/2

x1 <- runif(n)
y1 <- x1*1/4 + rnorm(n,0,.2)

x <- c(x1, rnorm(n))
y <- c(y1, rnorm(n))

dd <- data.frame(x=x,y=y)


cor.function <- function(x,ind){cor.test(x[ind,1],x[ind,2],
                                         method = "spearman")$estimate}

cor.function(dd)

sim1 <- boot::boot(dd[1:50,], cor.function  ,R=100)

sim2 <- boot::boot(dd[51:100,], cor.function  ,R=100)


hist(sim1$t)
hist(sim2$t)

db <- (sim1$t0 - sim2$t0)^2
se <-(sd(sim1$t)^2)+(sd(sim2$t)^2)
td <- db/se
p.chi <- pchisq(td, df = 1, lower.tail = F)
p.nd <- pnorm(sqrt(td), lower.tail = F)

sim1$t0
sim2$t0
p.chi
p.nd


s2.b <- .7
s1.b <- .6
s2.sd <- .05
s1.sd <- .04


db <- (s2.b-s1.b)^2
sd <-(s2.sd^2)+(s1.sd^2)
td <- db/sd
p.chi<- pchisq(td, df = 1,lower.tail = F) # Chi-square distribution with 1df
p.norm<- (pnorm( sqrt(td) , lower.tail = F))*2 # Normal distribution 

p.chi
p.norm


n <- 200

MGSmatrix <- as.data.frame(matrix(rnorm(n*161),ncol=161))
names(MGSmatrix) <- paste0("HG3A.",1:161)

pheno <- data.table(SCAPISid = paste0("sample",1:n),
                    age = runif(n,55,64),
                    Sex = sample(c("male","female"),n,replace=T))

pheno <- cbind(pheno,MGSmatrix)
head(pheno)
dim(pheno)



# Results checking
library(data.table)
library(tidyverse)


setwd('/Users/gabba126/Documents/PhD_projects/1.Sleep/Results/')

cutlast <- function(char,n){
  l <- nchar(char)
  a <- l-n+1
  return(substr(char,a,l))
}

res.nobmi <- fread('cor_all.var_mgs.tsv')
res.nobmi %>% filter(q.value<.05) %>% group_by(exposure) %>% count()

res.bmi <- fread('cor.bmi_all.var_mgs.tsv')
res.bmi %>% filter(q.value<.05) %>% group_by(exposure) %>% count()

res.ext <- fread('cor2_all.var_mgs.tsv')
res.ext %>% filter(q.value<.05) %>% group_by(exposure) %>% count()
res.ext %>% filter(q.value<.05 & rho>0) %>% group_by(exposure) %>% count()
res.ext %>% filter(q.value<.05 & rho<0) %>% group_by(exposure) %>% count()

mgs.t90 <- res.ext[q.value<.05 & exposure=="t90",MGS]
mgs.odi <- res.ext[q.value<.05 & exposure=="odi",MGS]

mgs.t90.pos <- res.ext[q.value<.05 & rho>0 & exposure=="t90",MGS]
mgs.odi.pos <- res.ext[q.value<.05 & rho>0 & exposure=="odi",MGS]

mgs.t90.pos[mgs.t90.pos %in% mgs.odi.pos]

mgs.t90.neg <- res.ext[q.value<.05 & rho<0 & exposure=="t90",MGS]
mgs.odi.neg <- res.ext[q.value<.05 & rho<0 & exposure=="odi",MGS]

mgs.t90.neg[mgs.t90.neg %in% mgs.odi.neg]

res.ext %>% filter(exposure=="t90" & q.value<0.05)

library(rio)
taxonomy <- import('HG3.A.7_tax.xlsx')
names(taxonomy)[1] <- "mgs"

res.ext[, mgs :=cutlast(MGS,9)]

res.ext <- merge(res.ext, taxonomy, all.x = T, by="mgs")

res.ext %>% filter(MGS %in% mgs.t90.neg[mgs.t90.neg %in% mgs.odi.neg]) %>% select(MGS,order,family) %>% unique %>% group_by(order) %>% count()

res.ext %>% filter(exposure=="t90" & q.value<0.05 & rho<0) %>% select(MGS,rho,q.value,order,family)

res.ext %>% filter(exposure=="t90" & q.value<0.05 & rho<0) %>% select(MGS,rho,q.value) %>% arrange(q.value)

res.ext %>% filter(exposure %in% c("t90","odi") & q.value<.05) %>% summarize(unique(MGS))

# Imputation analyssi 
res.imp <- fread('cor_ahi_imput_mgs.tsv')
res.imp[,MGS := gsub("_",".",MGS)]
res.imp <- merge(res.imp, unique(res.ext[,.(MGS,mgs)]),by.x="MGS",by.y="mgs",all.x=T, all.y = F)

res.imp[q_value<.10,] %>% arrange(q_value)
                                                                             

# medication sensitiviy analysis 

res.med <- fread('cor.med_all.var_mgs.tsv')
res.med %>% filter(exposure =="t90" & MGS %in% mgs.t90) %>% select(MGS,rho,p.value,q.value)
res.med %>% filter(exposure =="odi" & MGS %in% mgs.odi) %>% select(MGS,rho,p.value,q.value)

res.atb <- fread('corsaatb_all.var_mgs.tsv')
res.atb %>% filter(exposure =="t90" & MGS %in% mgs.t90) %>% select(MGS,rho,p.value,q.value)
res.atb %>% filter(exposure =="odi" & MGS %in% mgs.odi) %>% select(MGS,rho,p.value,q.value)

res.lung <- fread('corsalung_all.var_mgs.tsv')
res.lung %>% filter(exposure =="t90" & MGS %in% mgs.t90) %>% select(MGS,rho,p.value,q.value)
res.lung %>% filter(exposure =="odi" & MGS %in% mgs.odi) %>% select(MGS,rho,p.value,q.value)

# HB stratified 
res.hb <- fread('cor_hb_stratified.tsv')



# Results enrichment for GMM ####

res.ea.gmm <- fread('ea_GMM.tsv')

res.ea.gmm[q.value<.05,.(Name,exposure, direction, q.value, HL1)] %>% arrange(q.value)


# Results health outcomes 
res.ho <- fread('cor_mgs.gmm_bphb.tsv')

res.ho[q.value<0.05 & model =="OSA model",.(MGS_features, outcomes, rho, q.value)]

names.gmm <- unique(res.ea.gmm[,.(pathway,Name)])

res.ho.name <- merge(res.ho,names.gmm, by.x="MGS_features", by.y = "pathway", all.x=T, all.y=F)

res.ho.name[q.value<0.05 & model =="OSA model" & rho>0,.(MGS_features,Name,outcomes,rho,q.value),by=outcomes]

res.ho.name[q.value<0.05 & model =="OSA model" & rho<0,.(MGS_features,Name,outcomes,rho,q.value),by=outcomes]


res.ho.name[q.value<0.05 & model =="OSA and BMI adjusted" & rho>0,.(MGS_features,Name,outcomes,rho,q.value)]
res.ho.name[q.value<0.05 & model =="OSA and BMI adjusted" & rho<0,.(MGS_features,Name,outcomes,rho,q.value)]



res.ho[q.value<0.05 & model =="OSA model",.(MGS_features)]
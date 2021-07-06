# Script 7 - Enrichment analysis of modules among the MGS identified with model 1

# Gabriel Baldanzi 2021-07-05

# Loading packages 
pacman::p_load(data.table,ggplot2, tidyr, fgsea,rio)

# input and output folders 
input1 = "/home/baldanzi/Sleep_apnea/Results/"
input2 = "/home/baldanzi/Datasets/MGS/original/"
output = "/home/baldanzi/Sleep_apnea/Results/"
output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

# Importing results 
res.ahi <- fread(paste0(input1,"cor_ahi_mgs.tsv"))
res.bmi <- fread(paste0(input1,"cor_BMI_mgs.tsv"))
res.t90 <- fread(paste0(input1,"cor_t90_mgs.tsv"))

table.res = rbind(res.bmi,res.ahi,res.t90)

# List of modules ####
load(paste0(input2,'MGS_HG3A.keggModule2MGS.RData')) # object = MGS_HG3A.keggModule2MGS

#pathway enrichment analyses
MGS_HG3A.keggModule2MGS=MGS_HG3A.keggModule2MGS[lapply(MGS_HG3A.keggModule2MGS,length)>0] ## you can use sapply,rapply
## Compute maximum length
max.length <- max(sapply(MGS_HG3A.keggModule2MGS, length))
## Add NA values to list elements
MGS_HG3A.keggModule2MGS <- lapply(MGS_HG3A.keggModule2MGS, function(v) { c(v, rep(NA, max.length-length(v)))})
## Rbind
MGS_HG3A.keggModule2MGS=data.frame(do.call(rbind, MGS_HG3A.keggModule2MGS),stringsAsFactors = F)

modules=import(paste0(input2,"upugut03.keggModuleComp.percent.tsv")) # modules number and description 
modules=modules[,1:4]

fgsea.module.fun=function(res,outcome)
{
  RANK=NULL
  for(x in outcome) {
    
    print(x)
    xxx=as.data.frame(res[res$exposure%in%x,])
    
    lev=rownames(MGS_HG3A.keggModule2MGS)
    mylist=NULL
    for(i in lev)
    {
      eval(parse(text=paste("mylist$",i,"=xxx[xxx[,'mgs']%in%MGS_HG3A.keggModule2MGS[i,],'MGS']",sep="")))
    }
    
    cor=xxx$cor.coefficient
    names(cor)=xxx$MGS
    
    set.seed(123)
    xRANK <- fgsea(pathways = mylist, stats = rank(cor), eps=0,scoreType="pos", maxSize=1500, nPermSimple = 100000)
    xRANK=cbind(x,xRANK)
    RANK=rbind(RANK,xRANK)
  }
  
  RANK=data.frame(RANK,stringsAsFactors = F)  
  RANK$q.value=RANK$padj
  names(RANK)[1]="exposure"
  
  RANK=merge(RANK,modules, by.x="pathway",by.y="Module")
  
  return(RANK)
}

expos = c("BMI","ahi","t90")

res=fgsea.module.fun(table.res,expos)
fwrite(res,paste0(output,"ea_modules_rho.tsv"),sep = "\t")



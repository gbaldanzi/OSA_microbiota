# Annotating Results 
fdr.ahi <- res.ahi[q.value<0.05,][["MGS"]]
fdr.t90 <- res.t90[q.value<0.05,][["MGS"]]
fdr.bmi <- res.bmi[q.value<0.05,][["MGS"]]


# Create key 

splitted <- strsplit(dat2[,1],"[|]")

lastname <- lapply(splitted, function(x){
  a <- x[[length(x)]]
})

HG3A <- lapply(lastname,function(x){
  substr(x,4,nchar(x))
})

alltaxa <- data.table(fullname = dat2[,1], 
                      lastname = unlist(lastname),
                      HG3A = unlist(HG3A))


# Create annotation for Lefse plot 

lefse_lists = data.frame(node=c(alltaxa[HG3A %in% fdr.ahi,lastname],
                                alltaxa[HG3A %in% fdr.t90,lastname],
                                alltaxa[HG3A %in% fdr.bmi,lastname]),
                         color=c(rep('indianred3',length(fdr.ahi)),
                                 rep('blue',length(fdr.t90)),
                                 rep('violet',length(fdr.bmi))),
                         stringAsFactor=F)

library(dplyr)
library(ggtree)
library(ggplot2)

p2 <- clade.anno2(p, lefse_lists, alpha=0.3)


pdf("/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/lefseplot_step2.pdf")
p2
dev.off()


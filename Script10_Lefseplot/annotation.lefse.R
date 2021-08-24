# Annotating Results 
fdr.ahi <- res.ahi[q.value<0.05,][["MGS"]]
fdr.t90 <- res.t90[q.value<0.05,][["MGS"]]
fdr.bmi <- res.bmi[q.value<0.05,][["MGS"]]

fdr.all <- fdr.ahi[fdr.ahi %in% fdr.bmi & fdr.ahi %in% fdr.t90]

fdr.ahi.excl <- fdr.ahi[!fdr.ahi %in% fdr.bmi & !fdr.ahi %in% fdr.t90]
fdr.t90.excl <- fdr.t90[!fdr.t90 %in% fdr.bmi & !fdr.t90 %in% fdr.ahi]
fdr.bmi.excl <- fdr.bmi[!fdr.bmi %in% fdr.t90 & !fdr.bmi %in% fdr.ahi]

fdr.ahi.t90 <- fdr.ahi[fdr.ahi %in% fdr.t90 & !fdr.ahi %in% fdr.bmi]
fdr.ahi.bmi <- fdr.ahi[fdr.ahi %in% fdr.bmi & !fdr.ahi %in% fdr.t90]
fdr.t90.bmi <- fdr.t90[fdr.t90 %in% fdr.bmi & !fdr.t90 %in% fdr.ahi]

# Create key 

splitted <- strsplit(dat[,1],"[|]")

lastname <- lapply(splitted, function(x){
  a <- x[[length(x)]]
})

HG3A <- lapply(lastname,function(x){
  substr(x,4,nchar(x))
})

alltaxa <- data.table(fullname = dat[,1], 
                      lastname = unlist(lastname),
                      HG3A = unlist(HG3A))


# Create annotation for Lefse plot 

nodelist <- c(fdr.all , fdr.ahi.excl,fdr.t90.excl, fdr.bmi.excl,
              fdr.ahi.t90 , fdr.ahi.bmi , fdr.t90.bmi )

cls <- c(rep('gray80', nrow(alltaxa[HG3A %in% fdr.all,])), 
         rep('green', nrow(alltaxa[HG3A %in% fdr.ahi.excl,])), 
         rep('orange', nrow(alltaxa[HG3A %in% fdr.t90.excl,])), 
         rep('indianred3', nrow(alltaxa[HG3A %in% fdr.bmi.excl,])), 
         rep('darkgreen', nrow(alltaxa[HG3A %in% fdr.ahi.t90,])), 
         rep('darkblue', nrow(alltaxa[HG3A %in% fdr.ahi.bmi,])), 
         rep('darkred', nrow(alltaxa[HG3A %in% fdr.t90.bmi,])))

  lefse_lists = data.frame(node=alltaxa[HG3A %in% nodelist,lastname],
                         color=cls)

library(dplyr)
library(ggtree)
library(ggplot2)

source('clade.anno_edited.fun.R')

p2 <- clade.anno_edited(p, lefse_lists, alpha=0.5)


pdf("/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/lefseplot_step2.pdf")
p2
dev.off()


# Test - remove before saving final version #### 
test.list <- lefse_lists[1:4,]
test.list$stringAsFactor <- NULL
test.list$color <- c("indianred3", "indianred3", "orange", "green")
p3 = clade.anno_edited(p,test.list,alpha=.5, anno.x=NULL, anno.y=NULL)
pdf("/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/lefse.test.plot3.pdf")
p3
dev.off()




tree <- rtree(20)
p <- ggtree(tree, layout = 'circular')
p1 <- p + geom_hilight(node=12) + geom_hilight(node=6, fill="red")
p1

set_hilight_legend <- function(p, color, label, alpha=.5) {
  d <- data.frame(color=color, var=label, x=0, y=1, alpha=alpha)
   p + geom_rect(aes_(fill=~var, xmin=~x, xmax=~x, ymin=~y, ymax=~y), data=d, inherit.aes=F) +
       guides(fill=guide_legend(override.aes=list(fill=alpha(d$color, d$alpha))))
}

  cls = c("indianred3","orange", "green")
  names = c("BMI", "AHI", "T90")
  set_hilight_legend(p= p1, color = cls, label = names) + theme(legend.position="right")
  

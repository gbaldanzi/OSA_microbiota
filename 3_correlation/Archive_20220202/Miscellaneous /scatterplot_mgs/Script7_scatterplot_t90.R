# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-08-12

# This script will produce scatter plots for MGSs correlated to 
# T90 in the fully adjusted model 

rm(list = ls())
pacman::p_load(data.table,ggplot2)

# input and output folders 
input1 = "/home/baldanzi/Sleep_apnea/Results/"

#output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"
output.plot = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

# Import results 
res.t90 <- fread(paste0(input1,"cor2_t90_mgs.tsv"))
res.bmi <- fread(paste0(input1,"cor2_BMI_mgs.tsv"))

res.t90[,q.value:=round(q.value,3)]
res.bmi[,q.value:=round(q.value,3)]

# Select the relevant MGSs 
mgs.bmi <- res.bmi[q.value<.05,]
mgs.t90 <- res.t90[q.value<.05 & !MGS %in% mgs.bmi$MGS,MGS]

# Scatter plot  ####

# Import full data 
valid.t90 <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/valid.t90_MGS.shannon_Upp.rds")

# Restrict to the relevant MGS
names.not.mgs <- names(valid.t90)[-(grep("___",names(valid.t90)))]
a <- c(names.not.mgs,mgs.t90)
valid.t90 <- valid.t90[,a,with=F]


# Scatter plot
plot_scatterplot = function (var.y,var.x,data) {
  ggplot(data, aes_string(x=var.x)) +
    geom_point(aes_string(y=var.y)) +
    ylab(var.y) + xlab(var.x)+
    ggtitle(var.y)
}

myplots <- lapply(mgs.t90, plot_scatterplot, data = valid.t90, var.x="t90")

pdf(paste0(output.plot,"plot_scatter_mgs_t90.pdf"))
myplots
dev.off()

# Ranked Scatter plot
plot_scatterplot = function (var.y,var.x,data) {
  ggplot(data, aes_string(x=var.x)) +
    geom_point(aes_string(y=var.y)) +
    ylab(paste("Ranked", var.y)) + xlab(paste("Ranked",var.x)) +
    ggtitle(var.y)
}

names.ranks = paste0("r.",mgs.t90)
valid.t90[,(names.ranks) := lapply(.SD,rank), .SDcols = mgs.t90]
valid.t90[,r.t90 := rank(t90)]
myplots <- lapply(names.ranks, plot_scatterplot, data = valid.t90, var.x="r.t90")


pdf(paste0(output.plot,"plot_scatter_ranked_mgs_t90.pdf"))
myplots
dev.off()

# Scatter plot - log+1

    # Transform rel. abundance to ln+1
    valid.t90[,(mgs.t90) := lapply(.SD,function(x){log(x+1)}), .SDcols = mgs.t90]


    plot_scatterplot = function (var.y,var.x,data) {
      ggplot(data, aes_string(x=var.x)) +
        geom_point(aes_string(y=var.y)) +
        ylab(paste0("log(",var.y,"+1)")) + xlab(var.x)+
        ggtitle(var.y)
      }

    myplots <- lapply(mgs.t90, plot_scatterplot, data = valid.t90, var.x="t90")

    pdf(paste0(output.plot,"plot_scatter_log_mgs_t90.pdf"))
    myplots
    dev.off()

# Scatter plot - CLR ####

clr.count <- fread('/home/baldanzi/Datasets/MGS/clean/MGS_clr.transformed_4839_upp.tsv')

a = c("SCAPISid",mgs.t90)
clr.count <- clr.count[,a,with=F]


dades <- valid.t90[,grep("____",names(valid.t90)):=NULL]

dades <- merge(dades, clr.count, by="SCAPISid", all.x=T, all.y=F)

plot_scatterplot = function (var.y,var.x,data) {
  ggplot(data, aes_string(x=var.x)) +
    geom_point(aes_string(y=var.y)) +
    ylab(paste0("clr(",var.y,")")) + xlab(var.x)+
    ggtitle(var.y)
}


myplots <- lapply(mgs.t90, plot_scatterplot, data = dades, var.x="t90")

pdf(paste0(output.plot,"plot_scatter_clr_mgs_t90.pdf"))
myplots
dev.off()

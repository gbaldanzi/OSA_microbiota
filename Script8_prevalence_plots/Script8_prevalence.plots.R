# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi, 2021-09-06

# Last update: 2022-01-03

# This script will create plots to investigate non-linear 
# associations between MGS and AHI, T90, and ODI


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

  clean.y.axis <- function(y){
  #y <- gsub("Absicoccus_", "A. ", y)
  yy <- unlist(strsplit(y, split = "____"))
  y <- paste0(yy[2],"\n",yy[1])
  y <- gsub("Bifidobacterium_longum_subsp._longum", "B. longum subsp.\nlongum", y)
  y <- gsub("Blautia_","B. ",y)
  y <- gsub("Pediococcus_", "P. ", y)
  y <- gsub("_AM42_11","", y)
  y <- gsub("Staphylococcus_","S. ", y)
  y <- gsub("_sp"," sp", y)
  y <- gsub("_v"," v",y)
  y <- gsub("_AF"," AF",y)
  y <- gsub("_8L"," 8L",y)
  y <- gsub("_ur"," ur",y)
  y <- gsub("Gordonibacter","G. ",y)
  return(y)
  }
  
  cutlast <- function(char,n=9){
    l <- nchar(char)
    a <- l-n+1
    return(substr(char,a,l))
  }
  
  # Line plots of prevalence 
  pa <- paste0("pa_",mgs.fdr)
  pheno[,(pa):= lapply(.SD, decostand, method="pa"), .SDcols = mgs.fdr]
  
  # Heatmap ####
  res.matrix <- res.m2[MGS %in% mgs.fdr,]
  
  a = c("MGS", "exposure", "cor.coefficient")
  clust.data <- res.matrix[,a,with=F] %>% pivot_wider(id_cols = MGS, 
                                           names_from = exposure, 
                                           values_from= cor.coefficient)
  
  temp <- scale(clust.data[,c("ahi", "odi", "t90")])
  ord <- hclust(dist(temp, method="euclidean"), method="ward.D")$order
  ord <- rev(ord)
  
  res.matrix[,exposure:=toupper(exposure)]
  res.matrix[,exposure:= factor(exposure, levels = c("BMI","T90", "ODI", "AHI"))]
  res.matrix[q.value<.05 , sig := "+" ]
  res.matrix[,mgs := factor(mgs, levels = cutlast(clust.data$MGS[ord],9))]
  
 
  
  
  hm <- ggplot(res.matrix, aes(x=mgs , y= exposure, fill = cor.coefficient)) + 
    geom_tile(colour="white") + geom_text(aes(label=sig), size=1) + 
    #guides(fill=guide_legend(title="\u03C1")) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         guide=guide_legend(keywidth = .3, keyheight = .5)) +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(size=6), 
          axis.text.x = element_text(size=4.5, angle=45, hjust = 0), 
          axis.ticks = element_blank(), 
          legend.position = "right", 
          legend.text = element_text(size = 5),
          legend.title = element_blank(), 
          plot.margin = unit(c(1,0,7,.3), "mm")) + 
    scale_x_discrete(position = "top", expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0))

  #grobs <- ggplotGrob(hm)$grobs
  #legend <- get_legend(hm)
  
  #hm <- plot_grid(hm+theme(legend.position = "none"), legend, ncol = 2, 
 #                 rel_heights = c(1, 1.2), rel_widths = c(1,.06), axis="t")
  
  # PDF file to check if the Heatmap was constructed properly 
  ggsave(file="check_hm.pdf", plot=hm)
    
    
    # Prevalence by groups 

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
    
    line.plot.fun <- function(y.axis,dd=dades) {
      require(data.table)
      require(ggplot2)
      require(dplyr)
      require(Hmisc)
    
    r <- round(res.m2[exposure=="ahi" & MGS==y.axis,rho],3)
    
    p1 <-  ggplot(data=dd, aes_string(x="Group", y=y.axis, group="exposure")) + 
      geom_line(aes(colour=exposure), linetype="dashed", size=0.1) +
      geom_point(aes(colour=exposure),size=1.1,stroke=0, shape=15)+
      ggtitle(ggtitle(gsub("____","\n",clean.y.axis(y.axis)))) +
      #scale_x_discrete(labels=names(n.gr)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = .5, face = 'bold',size=5.5),
            #plot.subtitle = element_text(hjust = .5, face = 'bold', size=8),
            axis.text = element_text(size=5.3), 
            #axis.title.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=4.5, angle=45),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(r=5.5, l=3.5,t=2.5, b=3.5,unit="pt")) + 
            #panel.border = element_rect(linetype = "solid", colour = "black", size=0.5)) +
      scale_color_manual(values=c("orange","grey75","cornflowerblue"))
    
    return(p1)
    }
    
    # Creating the plots 
    if(dir.exists("lineplots")==F){dir.create("lineplots")}
    
    mgs.fdr <- clust.data$MGS[ord]
    
        lineplots_3x <- lapply(mgs.fdr, line.plot.fun, dd=dades)
    names(lineplots_3x) <- mgs.fdr
    
    lineplots_3x[[1]] <-  lineplots_3x[[1]] + 
      guides(color=guide_legend(title="groups by:")) +
      theme(legend.position = "right", 
            legend.text = element_text(size=5),
            legend.key.size = unit(0.5,"mm"), 
            legend.title = element_text(size=5), 
            legend.margin = margin(l=-6,r=0,t=-8,unit = "mm"),
            plot.margin = margin(r=5.5, l=2.5,t=3.5, b=4.5,unit="pt"))
                                                
    
    saveRDS(lineplots_3x, file = paste0(output.plot,'/lineplots/lineplots_3x.rds'))
    
    lineplots.merged <- plot_grid(plotlist=lineplots_3x, ncol=6, labels = NULL)
    
    final.plot <- plot_grid(lineplots.merged, hm, labels=NULL, ncol=1, rel_heights = c(8,1))
    
    
    
    ggsave("lineplot_hm_merged.pdf", plot = final.plot, path="lineplots/")
    ggsave("lineplots_merged.pdf", plot = lineplots.merged, path="lineplots/")
    
    pdf(file = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/lineplots_pdf.pdf", 
        width = 7, height = 10)
    final.plot
    dev.off()
    
    
    
   
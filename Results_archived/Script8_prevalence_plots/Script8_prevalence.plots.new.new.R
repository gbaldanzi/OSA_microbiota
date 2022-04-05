# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-15

# This script will create plots to investigate non-linear 
# associations between top 10 species and T90 and ODI

rm(list = ls())
pacman::p_load(data.table,tidyverse,vegan, cowplot,ggpubr)


  # Folders 
  input <-  "/home/baldanzi/Datasets/sleep_SCAPIS/"
  results.folder <-  "/home/baldanzi/Sleep_apnea/Results/"
  outupt.plot <- "/home/baldanzi/Sleep_apnea/Results/Plots"
  wrf <-  "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/"

  # Importing data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  

  # Import results to create plots 
  res.fm <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
  
  # Selecting the species associated with a q-value<.05 
  
  fdr.t90 <- res.fm[exposure=="t90" & q.value < .05, MGS]
  fdr.odi <- res.fm[exposure=="odi" & q.value < .05, MGS]
  fdr.both <- unique(c(fdr.t90,fdr.odi))
  
  ## Calculating prevalence and selecting the most prevalent species associated with T90 and ODI
  
  # T90
  prev <- paste0("prev_",fdr.t90)
  pheno[,(prev):= lapply(.SD, decostand, method="pa"), .SDcols = fdr.t90]
  prev.t90 <- lapply(prev,function(x){round((sum(pheno[[x]])/nrow(pheno))*100,1)})
  prev.t90 <- data.table(MGS=substr(prev,6,100) , prevalence = do.call(rbind,prev.t90))
  prev.t90.mgs <- prev.t90[order(-prevalence.V1),][1:10,]
  
  # ODI
  prev <- paste0("prev_",fdr.odi)
  pheno[,(prev):= lapply(.SD, decostand, method="pa"), .SDcols = fdr.odi]
  prev.odi <- lapply(prev,function(x){round((sum(pheno[[x]])/nrow(pheno))*100,1)})
  prev.odi <- data.table(MGS=substr(prev,6,100) , prevalence = do.call(rbind,prev.odi))
  prev.odi.mgs <- prev.odi[order(-prevalence.V1),][1:10,]
  
  # Box plot of abundance 
  box.plot.fun <- function(x.axis,y.axis) {
    ggplot(dades, aes_string(x=x.axis,y=y.axis)) +
      geom_line(size=0.2) +
      geom_point(size=2,stroke=0, shape=15) 
  }
  
  dades <- copy(pheno)
  dades <- dades[valid.t90 =="yes",]
  
  dades <- dades %>% group_by(odicat) %>% summarise_at(prev.odi.mgs$MGS,median)
  
  plot.list <- lapply(prev.odi.mgs$MGS, box.plot.fun, x.axis="odicat")
  
  p1 <- ggarrange(plotlist = plot.list, nrow=5,ncol=2)
  

     # Line plots of abundance ####
    
    # T90
  
    dades <- copy(pheno)
    dades <- dades[valid.t90 =="yes",]
    
    dades <- dades %>% group_by(t90cat) %>% summarise_at(prev.t90.mgs$MGS,median)
    dades[,1] <- c("T90=0","t1","t2","t3")
    names(dades)[1] <- "Group"
    
    m.dades <- as.matrix(dades[,-1])
    rownames(m.dades) <- dades$Group
    
    dades <- data.table(Group = row.names(m.dades) , 
                            apply(m.dades,2,function(x) log2(x)-log2(x[1])))
    
    dades <- pivot_longer(dades, cols = names(dades[,-1]),
                              names_to = "MGS", values_to = "abundance")
    
    dades$MGS <- gsub("____"," (",dades$MGS) 
    dades$MGS <- gsub("prev_","",dades$MGS)
    dades$MGS <- gsub("_"," ",dades$MGS)
    dades$MGS <- paste0(dades$MGS,")")
    dades$Group <- factor(dades$Group, levels = c("T90=0","t1","t2","t3"))
    
    colrs <- c("#E31A1C", "#BEAED4", "#E6AB02", "#FFFF99", "#386CB0", "#F0027F", 
                 "#BF5B17", "#666666", "#66A61E", "#A6CEE3")
    
    p1 <-  ggplot(data=dades, aes_string(x="Group", y="abundance", group="MGS")) + 
      geom_line(aes(colour=MGS), size=0.2) +
      geom_point(aes(colour=MGS),size=2,stroke=0, shape=15) +
      xlab("T90 groups") + ylab("log2 fold-change of relative abundance") +
      scale_color_manual(values=colrs) + 
      scale_x_discrete(expand=c(0.05,0.05)) +
      theme(legend.title = element_blank(), 
            legend.text = element_text(size=7),
            axis.title.y = element_text(size=9), 
            axis.text.y = element_text(size=8))
    
    ggsave(plot = p1, path=wrf, file="test.png")
    

    # ODI
    dades <- copy(pheno)
    dades <- dades[valid.t90 =="yes",]
    
    dades <- dades %>% group_by(odicat) %>% summarise_at(prev.odi.mgs$MGS,median)
    dades[,1] <- c("q1","q2","q3","q4")
    names(dades)[1] <- "Group"
    
    m.dades <- as.matrix(dades[,-1])
    rownames(m.dades) <- dades$Group
    
    dades <- data.table(Group = row.names(m.dades) , 
                        apply(m.dades,2,function(x) log2(x)-log2(x[1])))
    
    dades <- pivot_longer(dades, cols = names(dades[,-1]),
                          names_to = "MGS", values_to = "abundance")
    
    dades$MGS <- gsub("____"," (",dades$MGS) 
    dades$MGS <- gsub("prev_","",dades$MGS)
    dades$MGS <- gsub("_"," ",dades$MGS)
    dades$MGS <- paste0(dades$MGS,")")
    
    colrs <- c("#E31A1C", "#BEAED4", "#E6AB02", "#FFFF99", "#386CB0", "#F0027F", 
               "#BF5B17", "#666666", "#66A61E", "#A6CEE3")
    
    p2 <-  ggplot(data=dades, aes_string(x="Group", y="abundance", group="MGS")) + 
      geom_line(aes(colour=MGS), linetype="dashed", size=0.2) +
      geom_point(aes(colour=MGS),size=1.5,stroke=0, shape=15) +
      xlab("ODI groups") + ylab("log2 fold-change of relative abundance") +
      scale_color_manual(values=colrs) + 
      scale_x_discrete(expand=c(0.05,0.05)) +
      theme(legend.title = element_blank(), 
            legend.text = element_text(size=7),
            axis.title.y = element_text(size=9), 
            axis.text.y = element_text(size=8))

    ggsave(plot = p2, path=wrf, file="odi.png")

    
    pdf(file=paste0(wrf,"mer.pre.pdf"))
    plot_grid(p1,p2,labels = c("A","B"), label_size = 12, ncol=1)
    dev.off()
    
    library(RColorBrewer)
    n <- 10
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    pie(rep(1,10), col=colrs)
    pie(rep(1,30), col=col_vector[1:30])
    
    
   
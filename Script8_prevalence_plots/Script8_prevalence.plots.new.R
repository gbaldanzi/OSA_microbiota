# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last update: 2022-02-15

# This script will create plots to investigate non-linear 
# associations between top 10 species and T90 and ODI

rm(list = ls())
pacman::p_load(data.table,tidyverse,vegan, cowplot, RColorBrewer)


  # Folders 
  input <-  "/home/baldanzi/Datasets/sleep_SCAPIS/"
  results.folder <-  "/home/baldanzi/Sleep_apnea/Results/"
  outupt.plot <- "/home/baldanzi/Sleep_apnea/Results/Plots"
  wrf <-  "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/"

  # Importing data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  

  # Import results to create plots 
  res.fm <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
  res.fm[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  top.odi <- res.fm[exposure=="odi",]
  top.odi <- top.odi[,rho_abs:= abs(rho)]
  top.odi <- top.odi[order(-rho_abs),][1:10, MGS]
  
  top.t90 <- res.fm[exposure=="t90",]
  top.t90 <- top.t90[,rho_abs:= abs(rho)]
  top.t90 <- top.t90[order(-rho_abs),][1:10, MGS]
  
  top.mgs <- unique(c(top.odi,top.t90))
  
   # Line plots of prevalence 
    
    # T90
  
    dades <- copy(pheno)
    dades <- dades[valid.t90 =="yes",]
    
    pa <- paste0("pa_",top.t90)
    dades[,(pa):= lapply(.SD, decostand, method="pa"), .SDcols = top.t90]
    
    
    dades <- dades %>% group_by(t90cat) %>% summarise_at(pa,function(x){round((sum(x)/length(x))*100,1)})
    dades[,1] <- c("t0","t1","t2","t3")
    names(dades)[1] <- "Group"
    
    m.dades <- as.matrix(dades[,-1])
    rownames(m.dades) <- dades$Group
    
    dades <- data.table(Group = row.names(m.dades) , 
                            apply(m.dades,2,function(x) log2(x)-log2(x[1])))
    
    dades <- pivot_longer(dades, cols = names(dades[,-1]),
                              names_to = "MGS", values_to = "prevalence")
    
    dades$MGS <- gsub("____"," (",dades$MGS) 
    dades$MGS <- gsub("pa_","",dades$MGS)
    dades$MGS <- gsub("_"," ",dades$MGS)
    dades$MGS <- paste0(dades$MGS,")")
    
    colrs <- c("#E31A1C", "#BEAED4", "#E6AB02", "#FFFF99", "#386CB0", "#F0027F", 
                 "#BF5B17", "#666666", "#66A61E", "#A6CEE3")
    
    p1 <-  ggplot(data=dades, aes_string(x="Group", y="prevalence", group="MGS")) + 
      geom_line(aes(colour=MGS), linetype="dashed", size=0.2) +
      geom_point(aes(colour=MGS),size=1.5,stroke=0, shape=15) +
      xlab("T90 groups") + ylab("log2 fold-change") +
      scale_color_manual(values=colrs) + 
      scale_x_discrete(labels = c("t0"="T90=0"), expand=c(0.05,0.05)) +
      theme(legend.title = element_blank(), 
            legend.text = element_text(size=7),
            axis.title.y = element_text(size=9), 
            axis.text.y = element_text(size=8))
    
    #ggsave(plot = p1, path=wrf, file="test.png")
    
    # ODI
    dades <- copy(pheno)
    dades <- dades[valid.t90 =="yes",]
    
    pa <- paste0("pa_",top.odi)
    dades[,(pa):= lapply(.SD, decostand, method="pa"), .SDcols = top.odi]
    
    
    dades <- dades %>% group_by(odicat) %>% summarise_at(pa,function(x){round((sum(x)/length(x))*100,1)})
    dades[,1] <- c("q1","q2","q3","q4")
    names(dades)[1] <- "Group"
    
    m.dades <- as.matrix(dades[,-1])
    rownames(m.dades) <- dades$Group
    
    dades <- data.table(Group = row.names(m.dades) , 
                        apply(m.dades,2,function(x) log2(x)-log2(x[1])))
    
    dades <- pivot_longer(dades, cols = names(dades[,-1]),
                          names_to = "MGS", values_to = "prevalence")
    
    dades$MGS <- gsub("____"," (",dades$MGS) 
    dades$MGS <- gsub("pa_","",dades$MGS)
    dades$MGS <- gsub("_"," ",dades$MGS)
    dades$MGS <- paste0(dades$MGS,")")
    
    colrs <- c("#E31A1C", "#BEAED4", "#E6AB02", "#FFFF99", "#386CB0", "#F0027F", 
               "#BF5B17", "#666666", "#66A61E", "#A6CEE3")
    
    p2 <-  ggplot(data=dades, aes_string(x="Group", y="prevalence", group="MGS")) + 
      geom_line(aes(colour=MGS), linetype="dashed", size=0.2) +
      geom_point(aes(colour=MGS),size=1.5,stroke=0, shape=15) +
      xlab("ODI groups") + ylab("log2 fold-change") +
      scale_color_manual(values=colrs) + 
      scale_x_discrete(expand=c(0.05,0.05)) +
      theme(legend.title = element_blank(), 
            legend.text = element_text(size=7),
            axis.title.y = element_text(size=9), 
            axis.text.y = element_text(size=8))

    #ggsave(plot = p2, path=wrf, file="odi.png")

    
    pdf(file=paste0(wrf,"mer.pre.pdf"))
    plot_grid(p1,p2,labels = c("A","B"), label_size = 12, ncol=1)
    dev.off()
    
    library(RColorBrewer)
    n <- 10
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    pie(rep(1,10), col=colrs)
    pie(rep(1,30), col=col_vector[1:30])
    
    
   
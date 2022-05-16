# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script will create the Heatmap (Fig 4) containing the GMM enrichment
# analysis results and the associations with health outcomes 

# On the "x-axis", we displaysystolic blood pressure,
# diastolic blood pressure and Hb1Ac

# On the "y-axis", we display the gut metabolic modules (microbial pathways). 

rm(list=ls())


library(RColorBrewer)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)


  # Import the signature species 
  results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
  #results.folder <-'/Users/gabba126/Documents/PhD_projects/1.Sleep/Results/'
  #output.plot <- "/Users/gabba126/Documents/PhD_projects/1.Sleep/Manuscript/Figures/"
  wrf <- "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/"

# Function to fix strings
  
  fix.species.name.fun <- function(char){
    char <- gsub("AF22_8LB","AF22-8LB",char)
    char <- gsub("AM42_11","AM42-11",char)
    char <- gsub("TF06_15AC","TF06-15AC",char)
    char <- gsub("AF46_10NS","AF46-10NS",char)
    char <- gsub("AF36_15AT","AF36-15AT",char)
    char <- gsub("Marseille_P4005","Marseille-P4005", char)
    char <- gsub('[(]HG3A.0168[)]',"",char)
    char <- gsub("_"," ",char)
  }
      

  # GMM names 
   gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
   #gmm.names <- fread("/Users/gabba126/Documents/PhD_projects/1.Sleep/Manuscript/Figures/GMM_reference.csv")
  
   setnames(gmm.names,"Module","modules")
   gmm.names[,Name:=str_to_title(Name)]
   gmm.names[,Name:=gsub("Ii","II",Name)]
   gmm.names[,Name:=gsub("To ","to ",Name)]
   gmm.names[,Name:=gsub("-Off Phase","-off phase",Name)]

    
  # Import enrichment of GMM ####
    res.gmm = fread(paste0(results.folder,"ea_GMM.tsv"))
    res.gmm <- res.gmm[direction == "positive" & exposure == "t90",]
    
    res.gmm[,Name:=str_to_title(Name)]
    res.gmm[,Name:=gsub("Ii","II",Name)]
    res.gmm[,Name:=gsub("Iv","IV",Name)]
    res.gmm[,Name:=gsub("To ","to ",Name)]
    res.gmm[,Name:=gsub("Pay","pay",Name)]
    res.gmm[,Name:=gsub("Off Phase","off phase",Name)]
    
    res.gmm[q.value<.05, sig := "p.value.adj<.05"]
    res.gmm[q.value>.05, sig := "p.value.adj>.05"]
    
    res.gmm[, Name := factor(Name, levels = Name[order(NES)])]
    
    p1 <-  ggplot(res.gmm[pval<.20,], aes(x=NES,y=Name, fill = sig)) + 
      geom_bar(stat='identity') + ylab("Gut metabolic modules") + 
        labs(fill = bquote("p-value"[adj])) +
        scale_fill_manual(labels = c("<0.05", ">0.05"), values = c("darkred","gray80")) +
        theme_classic() + 
      theme(legend.title = element_text(size=8), 
            legend.text = element_text(size = 8), 
            axis.text.x = element_text(size=8), 
            axis.text.y = element_text(size=8.5, color = "black"))

  # Heatmap ####
    
  # Import results on the association with SBP/DBP/HbA1c 
    res.mgs.bp <- fread(paste0(results.folder, 'cor_mgs.gmm_bphb.tsv'))
    
    res.mgs.bp <- merge(res.mgs.bp, res.gmm[,.(pathway,Name)], by.x = "MGS_features", 
                        by.y = "pathway", all.x=T, all.y=F)
    
    res.mgs.bp[!is.na(Name), MGS_features := Name]
    
  

  # Association with SBP/DBP/HbA1c ####
    
    res.mgs.bp.ahi <- res.mgs.bp[model=="OSA model",]
    res.mgs.bp.bmi <- res.mgs.bp[model=="OSA and BMI adjusted",]
    
    res.mgs.bp.ahi[,outcomes := factor(outcomes, levels = c("SBP_Mean", "DBP_Mean", "Hba1cFormattedResult"),
                                       labels = c("SBP_AHI","DBP_AHI","HbA1c_AHI"))]
    
    res.mgs.bp.bmi[,outcomes := factor(outcomes, levels = c("SBP_Mean", "DBP_Mean", "Hba1cFormattedResult"),
                                       labels = c("SBP_BMI","DBP_BMI","HbA1c_BMI"))]
    
    
    res.mgs.bp <- rbind(res.mgs.bp.ahi, res.mgs.bp.bmi)
    
    #res.mgs.bp <- merge(res.mgs.bp, gmm.names[,.(modules,Name)], by.x = "MGS_features", by.y="modules", all.x=T)
    #res.mgs.bp[, MGS_features := Name]
    
    res.mgs.bp.cor <- res.mgs.bp[,.(MGS_features,outcomes,rho)] %>% 
      spread(key=outcomes, value = rho) 
    res.mgs.bp.q <- res.mgs.bp[,.(MGS_features,outcomes,q.value)] %>% 
      spread(key=outcomes, value = q.value)
    
    sig <- apply(res.mgs.bp.q,1,function(x) any(x<.05))
    
    res.mgs.bp.cor <- res.mgs.bp.cor[which(sig==T),]
    res.mgs.bp.q <- res.mgs.bp.q[which(sig==T),]
    
    setDF(res.mgs.bp.cor)
    setDF(res.mgs.bp.q)
    
    m <- res.mgs.bp.cor$MGS_features
    res.mgs.bp.cor$MGS_features <- NULL
    res.mgs.bp.cor <- as.matrix(res.mgs.bp.cor)
    rownames(res.mgs.bp.cor) <- m
    
    m <- res.mgs.bp.q$MGS_features
    res.mgs.bp.q$MGS_features <- NULL
    res.mgs.bp.q <- as.matrix(res.mgs.bp.q)
    rownames(res.mgs.bp.q) <- m
    
    
    # Associations with T90/ODI
    # Import results in the associations with T90/ODI
    
    res.bm <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
    res.bm[,exposure:=toupper(exposure)]
    
    res.bm.cor <- res.bm[exposure %in% c("T90","ODI") & MGS %in% m,.(MGS,exposure,rho)] %>%
      pivot_wider(id_cols=c(MGS), names_from = exposure, values_from = c(rho))
    
    res.bm.q <- res.bm[exposure %in% c("T90","ODI") & MGS %in% m,.(MGS,exposure,q.value)] %>%
      pivot_wider(id_cols=c(MGS), names_from = exposure, values_from = c(q.value))
    
    osa.gmm <- m[!m %in% res.bm.cor$MGS]
    mgs.pos <- res.bm.cor$MGS[res.bm.cor$T90>0]
    mgs.neg <- res.bm.cor$MGS[res.bm.cor$T90<0]
    
    order.features <- c(osa.gmm,mgs.pos,mgs.neg)
    
    # Create empty boxes for the osa.gmm in the res.bm.cor and res.bm.q matrix
    
    empty.box <- data.frame(MGS = osa.gmm, T90 = rep(0,length(osa.gmm)), ODI = rep(0,length(osa.gmm)))
    res.bm.cor <- rbind(empty.box,res.bm.cor)
    
    empty.box <- data.frame(MGS = osa.gmm, T90 = rep(1,length(osa.gmm)), ODI = rep(1,length(osa.gmm)))
    res.bm.q <- rbind(empty.box, res.bm.q)
    
    res.bm.cor <- res.bm.cor[match(order.features,res.bm.q$MGS ),]
    res.bm.q <- res.bm.q[match(order.features,res.bm.q$MGS ),]
    
    res.bm.cor$MGS <- res.bm.q$MGS <- NULL
    
    res.bm.cor <- as.matrix(res.bm.cor)
    res.bm.q <- as.matrix(res.bm.q)
    
    rownames(res.bm.cor) <- rownames(res.bm.q) <- order.features
    
    
    
  # Assert that the heatmap matrixes have the same order 
   assert.order <- function(x){
      return(x[match(order.features, rownames(x)),])
    }

   res.mgs.bp.cor <- assert.order(res.mgs.bp.cor)
   res.mgs.bp.q <- assert.order(res.mgs.bp.q)
  
  # Legend ####

  #col_fun3 = circlize::colorRamp2(c(-0.1, 0, 0.1), c(muted("darkcyan"),"white","goldenrod2"))
  

  #lgd3 = Legend(col_fun = col_fun3, at = c(-.1,-0.05,0,0.05,0.1), direction = "vertical", 
   #             title = expression(rho), legend_height = unit(25,"mm"), labels_gp = gpar(fontsize=7))
  
  # Annotations ####
  
  ha2 =  HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white", col = 0),
                                            labels = c("OSA adj.", "OSA+BMI adj."),
                                            labels_gp = gpar(col = "black", fontsize = 8)))
  
    
  # Fix names 
    noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(res.bm.cor)) ) )
    noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
    noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
    noms <- fix.species.name.fun(noms$V1) 
    rownames(res.bm.cor) <- noms
    
    noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(res.mgs.bp.cor)) ) )
    noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
    noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
    noms <- fix.species.name.fun(noms$V1) 
    rownames(res.mgs.bp.cor) <- noms
    
    rownames(res.bm.cor)[rownames(res.bm.cor) %in% osa.gmm] <- paste0(osa.gmm,"◊")
    rownames(res.mgs.bp.cor)[rownames(res.mgs.bp.cor) %in% osa.gmm] <- paste0(osa.gmm,"^1")
  

      # Create heatmap ####
  set.seed(7)
  ht_opt$TITLE_PADDING = unit(c(6, 6), "points")
  
  h1 <-  Heatmap(res.bm.cor, 
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if(res.bm.q[i, j] < 0.05) {
                     grid.text('*', x, y)
                   } },
                 column_names_rot =  40,
                 column_names_side = "bottom",
                 cluster_columns = T,
                 show_column_dend = F,
                 cluster_rows = F,
                 row_labels = c(expression("Lysine Degradation II"^1),
                                expression("Threonine Degradation I"^1),
                                rownames(res.bm.cor)[3:nrow(res.bm.cor)]),
                 
                 #column_labels = colnames(res.mgs.bp.cor),
                 column_names_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 8),
                 row_names_side = "left",
                 
                 row_split = factor(c(rep("gmm",length(osa.gmm)),
                                      rep("pos",length(mgs.pos)),
                                      rep("neg",length(mgs.neg))),
                                    levels = c("gmm","pos","neg")),
                 row_title = NULL,
                 
                 row_gap = unit(3,"mm"),
                 
                 
                 #column_gap = unit(3, "mm"),
                 #cluster_column_slices = F,
                 
                 #column_labels = rep(c("SBP", "DBP", "Hb1Ac"),2),
                 
                 #top_annotation = ha2,
                 #left_annotation = bar,
                 
                 col = circlize::colorRamp2( c(-0.1, 0, 0.1), c(muted("darkcyan") ,"white","goldenrod2")),
                 
                 heatmap_legend_param = list( title = expression(rho),
                                             labels_gp = gpar(fontsize=7),
                                            title_gp = gpar(fontsize=10),
                                           grid_width = unit(3, "mm"), 
                                          direction = "vertical"),
                 show_row_dend = F, 
                 width = unit(5*ncol(res.bm.cor),"mm")) +
  
   Heatmap(res.mgs.bp.cor, 
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(res.mgs.bp.q[i, j] < 0.05) {
                grid.text('*', x, y)
              } },
            column_names_rot =  40,
            column_names_side = "bottom",
            cluster_columns = T,
            show_column_dend = F,
            
            #column_labels = colnames(res.mgs.bp.cor),
            column_names_gp = gpar(fontsize = 8),
            row_names_gp = gpar(fontsize = 8),
            row_names_side = "left",
            
            column_title = "Health outcomes",
            column_title_gp = gpar(fontsize = 10,fontface='bold'),
            column_split = factor(c(rep("A", 3),
                                    rep("B",3)), 
                                  levels = c("A","B")),
            column_gap = unit(3, "mm"),
            cluster_column_slices = F,
            
            column_labels = rep(c("SBP", "DBP", "Hb1Ac"),2),
            
            top_annotation = ha2,
            #left_annotation = bar,
            
            col = circlize::colorRamp2( c(-0.1, 0, 0.1), c(muted("darkcyan") ,"white","goldenrod2")),
            
            #heatmap_legend_param = list( title = expression(rho),
             #                            labels_gp = gpar(fontsize=7),
              #                           title_gp = gpar(fontsize=10),
               #                          grid_width = unit(3, "mm"), 
                #                         direction = "vertical"),
            show_heatmap_legend = F,
            show_row_dend = F, 
            width = unit(5*ncol(res.mgs.bp.cor),"mm")) 

  

  p2 = grid.grabExpr(draw(h1))
  
  pfinal <- plot_grid(p1,p2, labels = c("a","b"))
  
  save_plot("Fig3.pdf", pfinal, base_height = 7.0, base_asp = 1.6)

  message("Saving plots")


 
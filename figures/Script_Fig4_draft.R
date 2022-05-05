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
library(tidyverse)
library(scales)


  # Import the signature species 
  #results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
  results.folder <-"/Users/gabba126/Documents/PhD_projects/1.Sleep/Manuscript/Figures/"
  output.plot <- "/home/baldanzi/Sleep_apnea/Results/Plots/"
  wrf <- "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/"

# Function to fix strings
  
    fix.subclass.name.fun <- function(char){
      char <- gsub("__PC_", " (PC)",char)
      char <- gsub("___","_",char)
      char <- gsub("__","_",char)
      char <- gsub("Drug", "Drug -", char)
      char <- gsub("_"," ",char)
      char <- gsub(" PI ", " (PI)",char)
      char <- gsub(" PE ", " (PE)",char)
      char <- gsub("Analgesics","Analgesics,", char)
      return(char)
    }
      

  # Import results on the enrichment of metabolites subclasses ####
   #ea.gmm.met <- fread(paste0(results.folder, "ea_gmm_subclass.tsv"))
  
  # GMM names 
   #gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
    gmm.names <- fread(paste0(results.folder,'GMM_reference.csv'))
  
   setnames(gmm.names,"Module","modules")
   gmm.names[,Name:=str_to_title(Name)]
   gmm.names[,Name:=gsub("Ii","II",Name)]
  
    
  # Import enrichment of GMM ####
    res.gmm = fread(paste0(results.folder,"ea_GMM.tsv"))
    res.gmm <- res.gmm[q.value<.05 & direction == "positive",]
    
    res.gmm[,Name:=str_to_title(Name)]
    res.gmm[,Name:=gsub("Ii","II",Name)]
    
    osa.gmm <- res.gmm[,pathway]
  
  # Import results on the association with SBP/DBP/HbA1c 
    res.mgs.bp <- fread(paste0(results.folder, 'cor_mgs.gmm_bphb.tsv'))
    
    res.mgs.bp <- res.mgs.bp[MGS_features %in% osa.gmm,]
  
  
  # Association with SBP/DBP/HbA1c ####
    
    res.mgs.bp.ahi <- res.mgs.bp[model=="OSA model",]
    res.mgs.bp.bmi <- res.mgs.bp[model=="OSA and BMI adjusted",]
    
    res.mgs.bp.ahi[,outcomes := factor(outcomes, levels = c("SBP_Mean", "DBP_Mean", "Hba1cFormattedResult"),
                                       labels = c("SBP_AHI","DBP_AHI","HbA1c_AHI"))]
    
    res.mgs.bp.bmi[,outcomes := factor(outcomes, levels = c("SBP_Mean", "DBP_Mean", "Hba1cFormattedResult"),
                                       labels = c("SBP_BMI","DBP_BMI","HbA1c_BMI"))]
    
    
    res.mgs.bp <- rbind(res.mgs.bp.ahi, res.mgs.bp.bmi)
    
    res.mgs.bp <- merge(res.mgs.bp, gmm.names[,.(modules,Name)], by.x = "MGS_features", by.y="modules", all.x=T)
    res.mgs.bp[, MGS_features := Name]
    
    res.mgs.bp.cor <- res.mgs.bp[,.(MGS_features,outcomes,rho)] %>% 
      spread(key=outcomes, value = rho) 
    res.mgs.bp.q <- res.mgs.bp[,.(MGS_features,outcomes,q.value)] %>% 
      spread(key=outcomes, value = q.value) 
    
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
  
    
  # Assert that the heatmap matrixes have the same order 
   # assert.order <- function(x){
  #    return(x[match(osa.gmm, rownames(x)),])
  #  }

#  res.mgs.bp.cor <- assert.order(res.mgs.bp.cor)
 # res.mgs.bp.q <- assert.order(res.mgs.bp.q)

  
  # Legend ####

  col_fun3 = circlize::colorRamp2(c(-0.1, 0, 0.1), c(muted("darkcyan"),"white","goldenrod2"))
  

  lgd3 = Legend(col_fun = col_fun3, at = c(-.1,-0.05,0,0.05,0.1), direction = "vertical", 
                title = expression(rho), legend_height = unit(25,"mm"), labels_gp = gpar(fontsize=7))
  
  # Annotations ####
  
  ha2 =  HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white", col = 0),
                                            labels = c("OSA adj.", "OSA+BMI adj."),
                                            labels_gp = gpar(col = "black", fontsize = 8)))
  
  # Bar plot of NES ####
  
  barplot.NES <- as.numeric(res.gmm[,NES])
  names(barplot.NES) <- res.gmm[,Name]
  
  bar = rowAnnotation("NES" = anno_barplot(barplot.NES,
                                                     gp = gpar(fill = "gray60", 
                                                               col="white", 
                                                               lwd = .05),
                                                     axis_param = list(gp=gpar(fontsize=7))),
                          annotation_name_side = "top",
                          annotation_name_gp = gpar(fontsize=9),
                          annotation_name_rot = 0,
                      empty = anno_empty(border = FALSE, width = unit(.5, "cm")))
                          #name = "NES (T90)")
                          #show_annotation_name = FALSE)
  

      # Create heatmap ####
  set.seed(7)
  ht_opt$TITLE_PADDING = unit(c(6, 6), "points")
  
  h1 <-  Heatmap(res.mgs.bp.cor, 
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
            left_annotation = bar,
            
            col = circlize::colorRamp2( c(-0.1, 0, 0.1), c(muted("darkcyan") ,"white","goldenrod2")),
            
            #heatmap_legend_param = list( title = expression(rho),
             #                            labels_gp = gpar(fontsize=7),
              #                           title_gp = gpar(fontsize=10),
               #                          grid_width = unit(3, "mm"), 
                #                         direction = "vertical"),
            show_row_dend = F, 
            width = unit(5*ncol(res.mgs.bp.cor),"mm")) 

  

  h1

  message("Saving plots")
  dev.off()
  
  # Draw ####
  #pdf(file = paste0(wrf, "Fig4_draft.pdf"), 
  pdf(file = paste0(results.folder,"Fig4_draft.pdf"), 
      width = 6, height = 3.2)
  set.seed(7)
  draw(h1, heatmap_legend_side = "right", show_heatmap_legend=F)
  draw(lgd3, x = unit(.98, "npc"), y = unit(.5, "npc"), just = c("right"))
  #draw(bar, x= unit(.371, "npc"), y= unit(.826, "npc"), width = unit(.163, "npc"), height = unit(.18, "npc"),
  #     angle=270, gp = gpar(lwd=.01, col="black"))
  dev.off()
  

 
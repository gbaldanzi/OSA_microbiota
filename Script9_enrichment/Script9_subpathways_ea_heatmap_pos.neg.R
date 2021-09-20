# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-16

# Heatmap on the Enrichment analysis of subpathways among the metabolites 
# correlated to the identified MGSs

rm(list=ls())


library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(tidyverse)

  # Important the relevant MGS for the heatmap 
  input = "/home/baldanzi/Sleep_apnea/Results/"

  res.m2 <- fread(paste0(input,"cor2_all.var_mgs.tsv"))

  res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]

  # Select the relevant MGSs 
  mgs.bmi <- res.m2[exposure =="BMI" & q.value<.05,MGS] # MGS correlated to BMI
  mgs.ahi <- res.m2[exposure =="ahi" & q.value<.05,MGS]
  mgs.t90 <- res.m2[exposure =="t90" & q.value<.05,MGS]

  # MGS correlated to either AHI or T90 but not to BMI (relevant mgs)
  mgs.rel <- res.m2 %>% filter(exposure =="ahi" | exposure == "t90") %>% 
    filter(q.value<.05) %>% filter(!MGS %in% mgs.bmi) %>% select(MGS)
  mgs.rel <- unique(mgs.rel$MGS)

  # Positive 

  # Import data
  res.table <- fread("/home/baldanzi/Sleep_apnea/Results/ea_subpathways_pos.tsv", data.table = F, sep = "\t")
  
  # Filter to only the pathways that are nominally associated with one MGS
  table.pathways <- res.table %>% group_by(pathway) %>% summarise(nr_p.05 = sum(pval<.05))
  hm.pathways <- table.pathways$pathway[table.pathways$nr_p.05>0]
  
  hm.matrix <- res.table %>% select(MGS, pathway, NES)%>% 
    filter(pathway %in% hm.pathways) %>% 
    spread(key=pathway, value = NES)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix <- as.matrix(hm.matrix)
  hm.matrix_pos <- t(hm.matrix)
  
  qvalues_pos <- res.table %>% select(MGS, pathway, padj)%>% 
    filter(pathway %in% hm.pathways) %>% 
    spread(key=pathway, value = padj)
  
  rownames(qvalues_pos) <- qvalues_pos$MGS
  qvalues_pos$MGS <- NULL
  
  qvalues_pos <- t(as.matrix(qvalues_pos))
  qvalues_pos[is.na(qvalues_pos)] <- 1
  
  
  # Negative 
  
  # Import data
  res.table <- fread("/home/baldanzi/Sleep_apnea/Results/ea_subpathways_neg.tsv", data.table = F, sep = "\t")
  
  # Filter to only the pathways that are nominally associated with one MGS
  table.pathways <- res.table %>% group_by(pathway) %>% summarise(nr_p.05 = sum(pval<.05))
  hm.pathways <- table.pathways$pathway[table.pathways$nr_p.05>0]
  
  hm.matrix <- res.table %>% select(MGS, pathway, NES)%>% 
    filter(pathway %in% hm.pathways) %>% 
    spread(key=pathway, value = NES)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix <- as.matrix(hm.matrix)
  hm.matrix_neg <- t(hm.matrix)
  
  qvalues_neg <- res.table %>% select(MGS, pathway, padj)%>% 
    filter(pathway %in% hm.pathways) %>% 
    spread(key=pathway, value = padj)
  
  rownames(qvalues_neg) <- qvalues_neg$MGS
  qvalues_neg$MGS <- NULL
  
  qvalues_neg <- t(as.matrix(qvalues_neg))
  qvalues_neg[is.na(qvalues_neg)] <- 1
  
  
  # annotation for AHI or T90 correlated MGS

      mgs.rel <- mgs.rel[mgs.rel %in% rownames(hm.matrix)]
      
      annotation <- data.table(mgs.rel = mgs.rel, 
                               correlated = as.character())
      
      annotation[mgs.rel %in% mgs.ahi, correlated := "AHI"]
      annotation[mgs.rel %in% mgs.t90, correlated := "T90"]
      annotation[mgs.rel %in% mgs.ahi & mgs.rel %in% mgs.t90, correlated := "Both"]
      
      annotation[,correlated:=factor(correlated, levels=c("AHI" , "Both" , "T90"))]
      
      setDF(annotation)
      
      annotation <- annotation[order(annotation$correlated),]
      
      rownames(annotation) <- annotation$mgs.rel 
      
      hm.matrix_neg <- hm.matrix_neg[,match(annotation$mgs.rel,colnames(hm.matrix_neg))]
      hm.matrix_pos <- hm.matrix_pos[,match(annotation$mgs.rel,colnames(hm.matrix_pos))]
      qvalues_pos <- qvalues_pos[,match(annotation$mgs.rel,colnames(qvalues_pos))]
      qvalues_neg <- qvalues_neg[,match(annotation$mgs.rel,colnames(qvalues_neg))]
      
      map.col = c("orange", "cornflowerblue", "red")
      names(map.col) = c("AHI", "T90", "Both")
      
      
      noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , colnames(hm.matrix_pos)) ) )
      noms <- paste0(noms$V1, " (", noms$V2, ")" )
      colnames(hm.matrix_pos) <- noms
      
      noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , colnames(hm.matrix_neg)) ) )
      noms <- paste0(noms$V1, " (", noms$V2, ")" )
      colnames(hm.matrix_neg) <- noms
      
      noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , rownames(annotation)) ) )
      noms <- paste0(noms$V1, " (", noms$V2, ")" )
      rownames(annotation) <- noms
      
  ha = HeatmapAnnotation(Correlation=annotation$correlated,which = "col",
                       name="Correlation",
                       col =  list(Correlation = c(map.col)),
                       annotation_legend_param = list(labels_gp = gpar(fontsize=4),
                                                      title_gp = gpar(fontsize=5),
                                                      grid_width = unit(2, "mm")),
                       annotation_label = "Correlation",
                       annotation_name_gp=gpar(fontsize = 5),
                       annotation_name_side = "right")
  
set.seed(1)
  h1 = Heatmap(hm.matrix_pos, 
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(qvalues_pos[i, j] < 0.05) {
                   grid.text('+', x, y)
                 } },
               column_names_rot =  90 , 
               column_names_side = "bottom",
               cluster_columns = F,
               show_column_dend = F,
               
               cluster_rows = TRUE,
               row_names_side = "left",
               show_row_dend = F,
               
               col = colorRamp2(c(0,1.5,3),c("white","indianred2", "red3")),
               name="NES (pos)",
               column_labels = colnames(hm.matrix_pos),
               column_names_gp = gpar(fontsize = 6),
               row_names_gp = gpar(fontsize = 6),
               #column_title_gp = gpar(fontsize = 7,fontface='bold'),
               heatmap_legend_param = list(labels_gp = gpar(fontsize=4),
                                           title_gp = gpar(fontsize=5),
                                           grid_width = unit(1.5, "mm")),
               top_annotation = ha) 
    
    h2 = Heatmap(hm.matrix_neg, 
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if(qvalues_neg[i, j] < 0.05) {
                     grid.text('+', x, y)
                   } },
            column_names_rot =  90,
            column_names_side = "bottom",
            cluster_columns = F,
            show_column_dend = F,
            
            cluster_rows = TRUE,
            row_names_side = "left",
            show_row_dend = F,
            
            col = colorRamp2(c(0,1.5,3),c("white","dodgerblue1","blue4" )),
            name="NES (neg)",
            column_labels = colnames(hm.matrix_neg),
            column_names_gp = gpar(fontsize = 6),
            row_names_gp = gpar(fontsize = 6),
            #column_title_gp = gpar(fontsize = 7,fontface='bold'),
            heatmap_legend_param = list(labels_gp = gpar(fontsize=4),
                                        title_gp = gpar(fontsize=5),
                                        grid_width = unit(1.5, "mm")),
            top_annotation = ha) 
    

#pdf(file = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/mgs_subpathway_ea_heatmap.pdf", width = 9, height = 5)
#set.seed(1)
#draw(h1, column_title = "Subpathways enrichment analysis",
#     column_title_gp = gpar(fontsize = 10, fontface="bold"),
#     heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend=T)

#dev.off()

png(filename =  "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/mgs_subpathway_ea_heatmap_pos.png", 
    width = 6, height = 8, units = 'in', res = 1500 )
set.seed(1)
draw(h1, column_title = "Subpathways enrichment analysis - positive correlations",
     column_title_gp = gpar(fontsize = 10, fontface="bold"),
     heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legend=T)

dev.off()

png(filename =  "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/mgs_subpathway_ea_heatmap_neg.png", 
    width = 6, height = 8, units = 'in', res = 1500 )
set.seed(1)
draw(h2, column_title = "Subpathways enrichment analysis - negative correlations",
     column_title_gp = gpar(fontsize = 10, fontface="bold"),
     heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legend=T)

dev.off()
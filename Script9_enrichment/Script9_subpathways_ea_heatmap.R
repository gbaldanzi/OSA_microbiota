# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-16

# Heatmap on the Enrichment analysis of subpathways among the metabolites 
# correlated to the identified MGSs

rm(list=ls())

# library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(tidyverse)

  # Import data
  res.table <- fread("/home/baldanzi/Sleep_apnea/Results/ea_subpathways.tsv", data.table = F)
  
  
  hm.matrix <- res.table %>% select(maintax, pathway, NES) %>% spread(key=pathway, value = NES)
  
  rownames(hm.matrix) <- hm.matrix$maintax
  hm.matrix$maintax <- NULL
  
  hm.matrix <- as.matrix(hm.matrix)
  
  # annotation for AHI or T90 correlated MGS
  input = "/home/baldanzi/Sleep_apnea/Results/"
  
    res.ahi <- fread(paste0(input,"cor2_ahi_mgs.tsv"))
    res.t90 <- fread(paste0(input,"cor2_t90_mgs.tsv"))
    res.bmi <- fread(paste0(input,"cor2_BMI_mgs.tsv"))
  
    res.ahi[,q.value:=round(q.value,3)]
    res.t90[,q.value:=round(q.value,3)]
    res.bmi[,q.value:=round(q.value,3)]
  
      # Select the relevant MGSs 
      mgs.bmi <- res.bmi[q.value<.05,]
      mgs.ahi <- res.ahi[q.value<.05 & !MGS %in% mgs.bmi$MGS,MGS]
      mgs.t90 <- res.t90[q.value<.05 & !MGS %in% mgs.bmi$MGS,MGS]
      
      mgs.rel <- unique(c(mgs.ahi,mgs.t90))
      
      #noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , mgs.rel) ) )
      #mgs.rel <- paste0(noms$V1, " (", noms$V2, ")" )
      
      mgs.rel <- mgs.rel[mgs.rel %in% rownames(hm.matrix)]
      
      annotation <- data.table(mgs.rel = mgs.rel, 
                               correlated = as.character())
      
      annotation[mgs.rel %in% mgs.ahi, correlated := "AHI"]
      annotation[mgs.rel %in% mgs.t90, correlated := "T90"]
      annotation[mgs.rel %in% mgs.ahi & mgs.rel %in% mgs.t90, correlated := "Both"]
  
      setDF(annotation)
      
      
      
      rownames(annotation) <- annotation$mgs.rel 
      
      map.col = c("#008000CC", "#466E82", "#FFA540")
      names(map.col) = c("AHI", "T90", "Both")
      
      
      noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , rownames(hm.matrix)) ) )
      noms <- paste0(noms$V1, " (", noms$V2, ")" )
      rownames(hm.matrix) <- noms
      
      noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , rownames(annotation)) ) )
      noms <- paste0(noms$V1, " (", noms$V2, ")" )
      rownames(annotation) <- noms
      
  ha = HeatmapAnnotation(Correlation=annotation$correlated,which = "row",
                       name="Correlation",
                       col = list(Correlation = c(map.col)),
                       annotation_legend_param = list(labels_gp = gpar(fontsize=4),
                                                      title_gp = gpar(fontsize=5),
                                                      grid_width = unit(2, "mm"),
                                                      nrow =1,
                                                      Correlation = list(nrow = 1)),
                       annotation_label = "Correlation",
                       annotation_name_gp=gpar(fontsize = 5),
                       annotation_name_rot = 90,
                       annotation_name_side = "bottom")



pdf(file = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/temp.heatmap.pdf", width = 9, height = 5)
set.seed(1)
h1 = Heatmap(hm.matrix, 
             column_names_rot =  90,
             column_names_side = "bottom",
             cluster_columns = TRUE,
             show_column_dend = FALSE,
             
             cluster_rows = TRUE,
             row_names_side = "right",
             row_dend_width = unit(1.5, "cm"),

             col = colorRamp2(c(0,2.3),c("gray95","darkred")),
             name="NES",
             column_labels = colnames(hm.matrix),
             column_names_gp = gpar(fontsize = 5),
             row_names_gp = gpar(fontsize = 5),
             #column_title_gp = gpar(fontsize = 7,fontface='bold'),
             heatmap_legend_param = list(labels_gp = gpar(fontsize=4),
                                         title_gp = gpar(fontsize=5),
                                         grid_width = unit(1.5, "mm"),
                                         direction = 'horizontal'),
             left_annotation = ha)


draw(h1, column_title = "Subpathways enrichment analysis",
     column_title_gp = gpar(fontsize = 10, fontface="bold"),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend=T)

dev.off()
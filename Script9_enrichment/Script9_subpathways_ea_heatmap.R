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

  # Import data
  res.table <- fread("/home/baldanzi/Sleep_apnea/Results/ea_subpathways.tsv", data.table = F, sep = "\t")
  
  # Filter to only the pathways that are nominaly associated with one MGS
  table.pathways <- res.table %>% group_by(pathway) %>% summarise(nr_p.05 = sum(pval<.05))
  hm.pathways <- table.pathways$pathway[table.pathways$nr_p.05>0]
  
  hm.matrix <- res.table %>% select(MGS, pathway, NES)%>% 
    filter(pathway %in% hm.pathways) %>% 
    spread(key=pathway, value = NES)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix <- as.matrix(hm.matrix)
  hm.matrix <- t(hm.matrix)
  
  # annotation for AHI or T90 correlated MGS
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
      
      #noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , mgs.rel) ) )
      #mgs.rel <- paste0(noms$V1, " (", noms$V2, ")" )
      
      mgs.rel <- mgs.rel[mgs.rel %in% colnames(hm.matrix)]
      
      annotation <- data.table(mgs.rel = mgs.rel, 
                               correlated = as.character())
      
      annotation[mgs.rel %in% mgs.ahi, correlated := "AHI"]
      annotation[mgs.rel %in% mgs.t90, correlated := "T90"]
      annotation[mgs.rel %in% mgs.ahi & mgs.rel %in% mgs.t90, correlated := "Both"]
  
      setDF(annotation)
      
      
      
      rownames(annotation) <- annotation$mgs.rel 
      
      map.col = c("#008000CC", "#466E82", "#FFA540")
      names(map.col) = c("AHI", "T90", "Both")
      
      
      noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , colnames(hm.matrix)) ) )
      noms <- paste0(noms$V1, " (", noms$V2, ")" )
      colnames(hm.matrix) <- noms
      
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
  h1 = Heatmap(hm.matrix, 
               column_names_rot =  90,
               column_names_side = "bottom",
               cluster_columns = TRUE,
               show_column_dend = T,
               
               cluster_rows = TRUE,
               row_names_side = "right",
               show_row_dend = F,
               
               col = colorRamp2(c(0,2.3),c("gray95","darkred")),
               name="NES",
               column_labels = colnames(hm.matrix),
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

png(filename =  "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/mgs_subpathway_ea_heatmap.png", 
    width = 6, height = 8, units = 'in', res = 1500 )
set.seed(1)
draw(h1, column_title = "Subpathways enrichment analysis",
     column_title_gp = gpar(fontsize = 10, fontface="bold"),
     heatmap_legend_side = "right", annotation_legend_side = "bottom", merge_legend=T)

dev.off()
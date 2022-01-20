# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2022-01-11

# Heatmap on the Enrichment analysis of GMM modules (Fig.5)

# Last update: 2022-01-11

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(tidyverse)
library(stringr)

  # input and output folders 
  input <-  "/home/baldanzi/Sleep_apnea/Results/"
  output <-  "/home/baldanzi/Sleep_apnea/Results/"

  #Importing
  res.pos <- fread(paste0(input,"ea_GMM_pos.tsv"))
  res.neg <- fread(paste0(input,"ea_GMM_neg.tsv"))

  # Import GMM names 
  gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  gmm.names[,Name:=str_to_title(Name)]
  gmm.names[,Name:=gsub("Ii","II",Name)]
  
  res.pos <- merge(res.pos,gmm.names,by.x="pathway",by.y="Module",all.x=T,all.y=F)
  res.neg <- merge(res.neg,gmm.names,by.x="pathway",by.y="Module",all.x=T,all.y=F)


  #Cleaning
  res.pos[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  res.neg[q.value>=0.001, q.value:=round(q.value, digits = 3)]
  
  res.pos[,exposure:=factor(exposure,levels = c("ahi","t90","odi","BMI"))]
  res.neg[,exposure:=factor(exposure,levels = c("ahi","t90","odi","BMI"))]
  
  #Matrix for heatmap 
  modules.hm.pos <- unique(res.pos[exposure %in% c("ahi","t90","odi") & q.value<.05,Name])
  modules.hm.neg <- unique(res.neg[exposure %in% c("ahi","t90","odi") & q.value<.05,Name])
  
  hm.pos <- res.pos %>% select(exposure, Name, NES)%>% 
    filter(Name %in% modules.hm.pos) %>% 
    spread(key=Name, value = NES)
  
  hm.neg <- res.neg %>% select(exposure, Name, NES)%>% 
    filter(Name %in% modules.hm.neg) %>% 
    spread(key=Name, value = NES)
  
  setDF(hm.pos)
  setDF(hm.neg)
  
  rownames(hm.pos) <- toupper(hm.pos$exposure)
  rownames(hm.neg) <- toupper(hm.neg$exposure)
  
  hm.pos$exposure <- hm.neg$exposure <- NULL
  
  hm.pos <- t(as.matrix(hm.pos))
  hm.neg <- t(as.matrix(hm.neg))
  
  
  # q-values 
  
  qvalues_pos <- res.pos %>% select(exposure, Name, q.value)%>% 
    filter(Name %in% modules.hm.pos) %>% 
    spread(key=Name, value = q.value)
  
  qvalues_neg <- res.neg %>% select(exposure, Name, q.value)%>% 
    filter(Name %in% modules.hm.neg) %>% 
    spread(key=Name, value = q.value)
  
  rownames(qvalues_pos) <- qvalues_pos$exposure
  rownames(qvalues_neg) <- qvalues_neg$exposure
  
  qvalues_pos$exposure <- qvalues_neg$exposure <-  NULL
  
  qvalues_pos <- t(as.matrix(qvalues_pos))
  qvalues_neg <- t(as.matrix(qvalues_neg))
  
  
  # Heatmap ####
  
  set.seed(1)
  h1 = Heatmap(hm.pos, 
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(qvalues_pos[i, j] < 0.05) {
                   grid.text('*', x, y,
                             gp = gpar(fontsize=6))
                 } },
               column_names_rot =  0 , 
               column_names_side = "top",
               column_names_centered = TRUE,
               cluster_columns = F,
               show_column_dend = F,
               
               #column_title_gp = gpar(fontsize=5),
               
               cluster_rows = TRUE,
               row_names_side = "left",
               show_row_dend = F,
               
               col = colorRamp2(c(.8,1.5,1.8),c("white","indianred3","red3"
                                             )),
               name="NES (pos)",
               #column_labels = colnames(hm.pos),
               column_names_gp = gpar(fontsize = 5),
               row_names_gp = gpar(fontsize = 5),
               #column_title_gp = gpar(fontsize = 7,fontface='bold'),
               heatmap_legend_param = list(labels_gp = gpar(fontsize=4),
                                           title_gp = gpar(fontsize=4),
                                           grid_width = unit(1, "mm"), 
                                           legend_height = unit(12.5, "mm"))) %v%
    
    Heatmap(hm.neg, 
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(qvalues_neg[i, j] < 0.05) {
                grid.text('*', x, y,
                          gp = gpar(fontsize=6))
              } },
            #column_names_rot =  0,
            #column_names_side = "top",
            show_column_names = FALSE,
            cluster_columns = F,
            show_column_dend = F,
            
            #column_title_gp = gpar(fontsize=5),
           
            cluster_rows = TRUE,
            row_names_side = "left",
            show_row_dend = F,
            
            col = colorRamp2(c(1.1,1.25,1.4),c("white","dodgerblue1",
                                            "blue4" 
                                          )),
            name="NES (neg)",
            #column_labels = colnames(hm.neg),
            #column_names_gp = gpar(fontsize = 4),
            row_names_gp = gpar(fontsize = 5),
            heatmap_legend_param = list(labels_gp = gpar(fontsize=4),
                                        title_gp = gpar(fontsize=4),
                                        legend_height = unit(12.5, "mm"),
                                        grid_width = unit(1, "mm"))) 
  
  
  
  # Drawing ####
  
  pdf(file = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/gmm_heatmap.pdf", 
      width = 4, height = 4)
  set.seed(1)
  draw(h1, heatmap_legend_side = "right", merge_legend=T)
  
  dev.off()
  
  
  
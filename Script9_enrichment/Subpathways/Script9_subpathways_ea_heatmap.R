# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script wil create the Heatmap on the Enrichment analysis of 
# subpathways among the metabolites 
# correlated to the signature species 

# The enrichment analysis results are imported from the Gutsy Atlas 
# reference: Dekkers, K. F. et al. An online atlas of human plasma metabolite 
# signatures of gut microbiome composition. 2021.12.23.21268179 
# https://www.medrxiv.org/content/10.1101/2021.12.23.21268179v1 (2021) 
# doi:10.1101/2021.12.23.21268179.

# Last update: 2022-02-03

rm(list=ls())


library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(tidyverse)


  # Import the signature species 
  results.folder <-  "/home/baldanzi/Sleep_apnea/Results/"
  output.plot <- "/home/baldanzi/Sleep_apnea/Results/Plots"

# Function to select the last characters of a string 
  cutlast <- function(char,n){
    l <- nchar(char)
    a <- l-n+1
    return(substr(char,a,l))
  }


  fix.species.name.fun <- function(char){
    char <- gsub("AF22_8LB","AF22-8LB",char)
    char <- gsub("AM42_11","AM42-11",char)
    char <- gsub("_"," ",char)
  }

  # Results from the MGS-AHI/T90/BMI correlation - Full model 

    res.m2 <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
    res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]
    
    Mgs.mgs <- unique(res.m2[,.(MGS,mgs)])

  # Signature species of OSA 

  mgs.m2 = readRDS(paste0(results.folder, "mgs.m2.rds"))
  mgs.rel <- unique(do.call('c',mgs.m2))
  mgs.rel.HG3A <- Mgs.mgs[MGS %in% mgs.rel, mgs]
  
  
  # Import enrichment analysis results GUTSY Atlas 
  ea_gutsy <- fread("/home/baldanzi/Datasets/gutsy_atlas/Supplementary_Table_6.tsv")
  
  names(ea_gutsy) <- gsub("-","_",names(ea_gutsy))
  names(ea_gutsy) <- gsub(" ","_",names(ea_gutsy))
  
  ea_gutsy[,mgs:=paste0("HG3A.",cutlast(Metagenomic_species,4))]
  
  ea_gutsy <- ea_gutsy[mgs %in% mgs.rel.HG3A,]
  ea_gutsy <- merge(ea_gutsy, Mgs.mgs, by ="mgs", all.x=T, all.y=F)
  
  
  # Increased and decreased species in OSA
  
  mgs.rel <- Mgs.mgs[MGS %in% mgs.rel & mgs %in% ea_gutsy$mgs,MGS]

  a <- c("ahi", "t90", "odi")

  mgs.increased <- unique(res.m2[cor.coefficient>0 & MGS %in% mgs.rel & exposure %in% a , MGS])
  mgs.decreased <- unique(res.m2[cor.coefficient<0 & MGS %in% mgs.rel & exposure %in% a , MGS])

  
# Enriched pathways/metabolite groups in the POSITIVE correlations ####

  res.table <- ea_gutsy[Direction=="Positive"]
  
  # Filter to pathways/metabolite groups that are FDR associated with at least one signature species 
  
  table.pathways <- res.table %>% group_by(Metabolite_subclass) %>% 
    summarise(nr_p.05 = sum(q_value<.05))
  hm.pathways <- table.pathways$Metabolite_subclass[table.pathways$nr_p.05>0]

  hm.pathways <- hm.pathways[-which(hm.pathways=="Food Component/Plant")]

  # Matrix results for the heatmap 
  hm.matrix <- res.table %>% select(MGS, Metabolite_subclass, Estimate)%>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = Estimate)
  
  setDF(hm.matrix)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix_pos <- t(as.matrix(hm.matrix))
  hm.matrix_pos[is.na(hm.matrix_pos)] <- 0

  rownames(hm.matrix_pos) <- gsub("__","_",rownames(hm.matrix_pos))
  rownames(hm.matrix_pos) <- gsub("_"," ",rownames(hm.matrix_pos))
  rownames(hm.matrix_pos) <- gsub("Metabolism","Metab.",rownames(hm.matrix_pos))

  # q-values to produce the "*" on the heatmap
  qvalues_pos <- res.table %>% select(MGS, Metabolite_subclass, q_value)%>% 
      filter(Metabolite_subclass %in% hm.pathways) %>% 
      spread(key=Metabolite_subclass, value = q_value)
  
  setDF(qvalues_pos)

  rownames(qvalues_pos) <- qvalues_pos$MGS
  qvalues_pos$MGS <- NULL

  qvalues_pos <- t(as.matrix(qvalues_pos))
  qvalues_pos[is.na(qvalues_pos)] <- 1


 # Enriched pathways/metabolite groups in the NEGATIVE correlations ####
 
  res.table <- ea_gutsy[Direction=="Negative"]

  # Filter to pathways/metabolite groups that are FDR associated with at least one signature species 
  table.pathways <- res.table %>% group_by(Metabolite_subclass) %>% 
    summarise(nr_p.05 = sum(q_value<.05))
  hm.pathways <- table.pathways$Metabolite_subclass[table.pathways$nr_p.05>0]
  
  # Matrix results for the heatmap 
    hm.matrix <- res.table %>% select(MGS, Metabolite_subclass, Estimate)%>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = Estimate)

  setDF(hm.matrix)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix_neg <- t(as.matrix(hm.matrix))
  hm.matrix_neg[is.na(hm.matrix_neg)] <- 0

  rownames(hm.matrix_neg) <- gsub("__","_",rownames(hm.matrix_neg))
  rownames(hm.matrix_neg) <- gsub("_"," ",rownames(hm.matrix_neg))
  rownames(hm.matrix_neg) <- gsub("Metabolism","Metab.",rownames(hm.matrix_neg))
  rownames(hm.matrix_neg) <- gsub(" (PC)","",rownames(hm.matrix_neg))

  # q-values to produce the "*" on the heatmap
    qvalues_neg <- res.table %>% select(MGS, Metabolite_subclass, q_value) %>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = q_value)
  
  setDF(qvalues_neg)

  rownames(qvalues_neg) <- qvalues_neg$MGS
  qvalues_neg$MGS <- NULL

  qvalues_neg <- t(as.matrix(qvalues_neg))
    qvalues_neg[is.na(qvalues_neg)] <- 1


# Annotation for increased or decreased species in OSA 

  mgs.rel <- mgs.rel[mgs.rel %in% rownames(hm.matrix)]

  annotation <- data.table(mgs.rel = mgs.rel, 
                         correlated = vector(mode="character", 
                                             length=length(mgs.rel)))

  annotation[mgs.rel %in% mgs.increased, correlated := "increased"]
  annotation[mgs.rel %in% mgs.decreased, correlated := "decreased"]
  
  annotation[,correlated:=factor(correlated, levels=c("increased" , "decreased"))]
  
setDF(annotation)

  annotation_increased <- annotation[mgs.rel %in% mgs.increased,]
  annotation_decreased <- annotation[mgs.rel %in% mgs.decreased,]

  annotation_increased <- annotation_increased[order(annotation_increased$correlated),]
  annotation_decreased <- annotation_decreased[order(annotation_decreased$correlated),]

  n.increased = nrow(annotation_increased)
  n.decreased = nrow(annotation_decreased)

  annotation <- rbind(annotation_increased,annotation_decreased)

  rownames(annotation) <- annotation$mgs.rel 

  # Assert that the heatmap matrix and the annotation vector have in the same order
  hm.matrix_neg <- hm.matrix_neg[,match(annotation$mgs.rel,colnames(hm.matrix_neg))]
  hm.matrix_pos <- hm.matrix_pos[,match(annotation$mgs.rel,colnames(hm.matrix_pos))]

  # Assert that the q-values matrix and the annotation vector have in the same order
  qvalues_pos <- qvalues_pos[,match(annotation$mgs.rel,colnames(qvalues_pos))]
  qvalues_neg <- qvalues_neg[,match(annotation$mgs.rel,colnames(qvalues_neg))]

  map.col = c("red", "blue")
  names(map.col) = c("increased", "decreased")

  # Fix species names to better visualization in the heatmap
  noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , colnames(hm.matrix_pos)) ) )
  noms <- paste0(noms$V1, " (", noms$V2, ")" )
  noms <- fix.species.name.fun(noms) 
  colnames(hm.matrix_pos) <- noms

  noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , colnames(hm.matrix_neg)) ) )
  noms <- paste0(noms$V1, " (", noms$V2, ")" )
  noms <- fix.species.name.fun(noms) 
  colnames(hm.matrix_neg) <- noms

  noms <- as.data.frame(do.call(rbind, strsplit( split = "____" , rownames(annotation)) ) )
  noms <- paste0(noms$V1, " (", noms$V2, ")" )
  noms <- fix.species.name.fun(noms) 
  rownames(annotation) <- noms

  # Create object with the annotations for the Heatmap 
  ha = HeatmapAnnotation(Correlation=annotation$correlated,which = "col",
                                 name="Phenotype",
                                show_legend = FALSE,
                                 simple_anno_size = unit(1, "mm"),
                                 col =  list(Correlation = c(map.col)),
                                 annotation_legend_param = list(labels_gp = gpar(fontsize=6),
                                                                title_gp = gpar(fontsize=7),
                                                                grid_width = unit(1, "mm")),
                                 annotation_label = "",
                                 annotation_name_gp=gpar(fontsize = 7),
                                 annotation_name_side = "left")

  # Create heatmap
  set.seed(1)
  h1 = Heatmap(hm.matrix_pos, 
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(qvalues_pos[i, j] < 0.05) {
                 grid.text('*', x, y)
               } },
             column_names_rot =  90 , 
             column_names_side = "bottom",
             cluster_columns = T,
             show_column_dend = F,
             
             column_split = factor(c(rep("Increased abundance with OSA",n.increased),
                                     rep("Decreased abundance with OSA",n.decreased))),
             column_title_gp = gpar(fontsize=9),
             column_gap = unit(2,"mm"),
             
             cluster_rows = TRUE,
             row_names_side = "left",
             show_row_dend = F,
             
             col = colorRamp2(c(0,1.5,3),c("white","indianred2", "red3")),
             name="NES (pos)",
             column_labels = colnames(hm.matrix_pos),
             column_names_gp = gpar(fontsize = 8),
             row_names_gp = gpar(fontsize = 8),
             #column_title_gp = gpar(fontsize = 7,fontface='bold'),
             heatmap_legend_param = list(labels_gp = gpar(fontsize=4),
                                         title_gp = gpar(fontsize=6),
                                         grid_width = unit(1.5, "mm")),
             top_annotation = ha) %v% # The symbol '%v%' means one heatmap over the other
  
  Heatmap(hm.matrix_neg, 
          cell_fun = function(j, i, x, y, w, h, fill) {
            if(qvalues_neg[i, j] < 0.05) {
              grid.text('*', x, y)
            } },
          column_names_rot =  40,
          column_names_side = "bottom",
          cluster_columns = T,
          show_column_dend = F,
          
          column_split = factor(c(rep("Increased abundance\nwith sleep apnea",n.increased),rep("Decreased abudance\nwith sleep apnea",n.decreased))),
          column_title_gp = gpar(fontsize=7),
          column_gap = unit(2,"mm"),
          
          cluster_rows = TRUE,
          row_names_side = "left",
          show_row_dend = F,
          
          col = colorRamp2(c(0,1.5,3),c("white","dodgerblue1","blue4" )),
          name="NES (neg)",
          column_labels = colnames(hm.matrix_neg),
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 8),
          #column_title_gp = gpar(fontsize = 7,fontface='bold'),
          heatmap_legend_param = list(labels_gp = gpar(fontsize=4),
                                      title_gp = gpar(fontsize=6),
                                      grid_width = unit(1.5, "mm"))) 




  # Save the plot in pdf 
  pdf(file = paste0(output.plot, "mgs_subpathway_ea_heatmap_gutsy.pdf"), 
    width = 10, height = 8)
  set.seed(1)
  draw(h1, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend=T)

  dev.off()

  # Save the plot in png 
  png(filename =  "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/mgs_subpathway_ea_heatmap_gutsy.png", 
  width = 9, height = 8, units = 'in', res = 1500 )
  set.seed(1)
  draw(h1, column_title = "Enrichment analysis for group of metabolites",
     column_title_gp = gpar(fontsize = 12, fontface="bold"),
     heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend=T)

  dev.off()


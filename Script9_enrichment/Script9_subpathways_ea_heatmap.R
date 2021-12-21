# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-08-16

# Heatmap on the Enrichment analysis of subpathways among the metabolites 
# correlated to the identified MGSs



# Last update: 2021-12-15

rm(list=ls())


library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(tidyverse)

# Import results from enrichment analysis (positive correlations)
# Import data
res.table <- fread("/home/baldanzi/Sleep_apnea/Results/ea_subpathways_pos.tsv", data.table = F, sep = "\t")

# Important the relevant MGS for the heatmap 
input = "/home/baldanzi/Sleep_apnea/Results/"

# results from the MGS-AHI/T90/BMI correlation - Full model 

res.m2 <- fread(paste0(input,"cor2_all.var_mgs.tsv"))
res.m2[q.value>=0.001, q.value:=round(q.value, digits = 3)]

# Slee apnea signature mgs 

mgs.m2 = readRDS(paste0(input, "mgs.m2.rds"))
mgs.rel <- unique(do.call('c',mgs.m2))
mgs.rel <- mgs.rel[mgs.rel %in% res.table$MGS]

#mgs.ahi <- mgs.m2$mgs.fdr.ahi
#mgs.t90 <- mgs.m2$mgs.fdr.t90
#mgs.odi <- mgs.m2$mgs.fdr.odi

#mgs.ahi.increased <- res.m2[cor.coefficient>0 & MGS %in% mgs.ahi & exposure=="ahi", MGS]
#mgs.ahi.decreased <- res.m2[cor.coefficient<0 & MGS %in% mgs.ahi & exposure=="ahi", MGS]

#exclusive.mgs.t90 <- mgs.t90[!mgs.t90 %in% mgs.ahi]

#mgs.t90.increased <- res.m2[cor.coefficient>0 & MGS %in% exclusive.mgs.t90 & exposure=="t90", MGS]
#mgs.t90.decreased <- res.m2[cor.coefficient<0 & MGS %in% exclusive.mgs.t90 & exposure=="t90", MGS]

#mgs.increased <- c(mgs.ahi.increased,mgs.t90.increased)
#mgs.decreased <- c(mgs.ahi.decreased, mgs.t90.decreased)

  a <- c("ahi", "t90", "odi")

  mgs.increased <- unique(res.m2[cor.coefficient>0 & MGS %in% mgs.rel & exposure %in% a , MGS])
  mgs.decreased <- unique(res.m2[cor.coefficient<0 & MGS %in% mgs.rel & exposure %in% a , MGS])

# Positive 


# Filter to only the pathways that are FDR associated with one MGS
table.pathways <- res.table %>% group_by(pathway) %>% summarise(nr_p.05 = sum(padj<.05))
hm.pathways <- table.pathways$pathway[table.pathways$nr_p.05>0]

  hm.pathways <- hm.pathways[-which(hm.pathways=="Partially_Characterized_Molecules")]
  hm.pathways <- hm.pathways[-which(hm.pathways=="Food_Component_Plant")]


hm.matrix <- res.table %>% select(MGS, pathway, NES)%>% 
  filter(pathway %in% hm.pathways) %>% 
  spread(key=pathway, value = NES)

rownames(hm.matrix) <- hm.matrix$MGS
hm.matrix$MGS <- NULL

hm.matrix <- as.matrix(hm.matrix)
hm.matrix_pos <- t(hm.matrix)
hm.matrix_pos[is.na(hm.matrix_pos)] <- 0

rownames(hm.matrix_pos) <- gsub("__","_",rownames(hm.matrix_pos))
rownames(hm.matrix_pos) <- gsub("_"," ",rownames(hm.matrix_pos))
rownames(hm.matrix_pos) <- gsub("Metabolism","Metab.",rownames(hm.matrix_pos))


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

# Filter to only the pathways that are FDR associated with one MGS
  table.pathways <- res.table %>% group_by(pathway) %>% summarise(nr_p.05 = sum(padj<.05))
  hm.pathways <- table.pathways$pathway[table.pathways$nr_p.05>0]
  
  
  hm.pathways <- hm.pathways[-which(hm.pathways=="Partially_Characterized_Molecules")]

  
  hm.matrix <- res.table %>% select(MGS, pathway, NES)%>% 
    filter(pathway %in% hm.pathways) %>% 
    spread(key=pathway, value = NES)

  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL

  hm.matrix <- as.matrix(hm.matrix)
  hm.matrix_neg <- t(hm.matrix)
  hm.matrix_neg[is.na(hm.matrix_neg)] <- 0

  rownames(hm.matrix_neg) <- gsub("__","_",rownames(hm.matrix_neg))
  rownames(hm.matrix_neg) <- gsub("_"," ",rownames(hm.matrix_neg))
  rownames(hm.matrix_neg) <- gsub("Metabolism","Metab.",rownames(hm.matrix_neg))
  rownames(hm.matrix_neg) <- gsub("Leucine","Leucine,",rownames(hm.matrix_neg))
  rownames(hm.matrix_neg) <- gsub("Glycolysis","Glycolysis,",rownames(hm.matrix_neg))

  rownames(hm.matrix_neg) <- gsub("PE","",rownames(hm.matrix_neg))
  rownames(hm.matrix_neg) <- gsub("PC","",rownames(hm.matrix_neg))


qvalues_neg <- res.table %>% select(MGS, pathway, padj)%>% 
  filter(pathway %in% hm.pathways) %>% 
  spread(key=pathway, value = padj)

rownames(qvalues_neg) <- qvalues_neg$MGS
qvalues_neg$MGS <- NULL

  qvalues_neg <- t(as.matrix(qvalues_neg))
    qvalues_neg[is.na(qvalues_neg)] <- 1


# annotation for increased or decreased MGS in sleep apnea 

  mgs.rel <- mgs.rel[mgs.rel %in% rownames(hm.matrix)]

  annotation <- data.table(mgs.rel = mgs.rel, 
                         correlated = vector(mode="character", 
                                             length=length(mgs.rel)))

  #annotation[mgs.rel %in% mgs.ahi, correlated := "AHI"]
  #annotation[mgs.rel %in% mgs.t90, correlated := "T90"]
  #annotation[mgs.rel %in% mgs.ahi & mgs.rel %in% mgs.t90, correlated := "Both"]

  #annotation[,correlated:=factor(correlated, levels=c("AHI" , "Both" , "T90"))]
  
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


  hm.matrix_neg <- hm.matrix_neg[,match(annotation$mgs.rel,colnames(hm.matrix_neg))]
  hm.matrix_pos <- hm.matrix_pos[,match(annotation$mgs.rel,colnames(hm.matrix_pos))]


qvalues_pos <- qvalues_pos[,match(annotation$mgs.rel,colnames(qvalues_pos))]
qvalues_neg <- qvalues_neg[,match(annotation$mgs.rel,colnames(qvalues_neg))]

  map.col = c("red", "blue")
  names(map.col) = c("increased", "decreased")


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
             
             column_split = factor(c(rep("Increased abundance\nwith sleep apnea",n.increased),rep("Decreased abundance\nwith sleep apnea",n.decreased))),
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
             top_annotation = ha) %v%
  
  Heatmap(hm.matrix_neg, 
          cell_fun = function(j, i, x, y, w, h, fill) {
            if(qvalues_neg[i, j] < 0.05) {
              grid.text('*', x, y)
            } },
          column_names_rot =  45,
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





pdf(file = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/mgs_subpathway_ea_heatmap.pdf", 
    width = 9, height = 8)
set.seed(1)
draw(h1, column_title = "Enrichment analysis for group of metabolites",
     column_title_gp = gpar(fontsize = 12, fontface="bold"),
     heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend=T)

dev.off()

png(filename =  "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/mgs_subpathway_ea_heatmap.png", 
width = 9, height = 8, units = 'in', res = 1500 )
set.seed(1)
draw(h1, column_title = "Enrichment analysis for group of metabolites",
     column_title_gp = gpar(fontsize = 12, fontface="bold"),
     heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend=T)

dev.off()


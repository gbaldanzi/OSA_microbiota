# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script wil create the Heatmap on the Enrichment analysis of 
# subpathways among the metabolites 
# correlated to the species associated with T90 and/or ODI in the full model 

# The enrichment analysis results are imported from the Gutsy Atlas 
# reference: Dekkers, K. F. et al. An online atlas of human plasma metabolite 
# signatures of gut microbiome composition. 2021.12.23.21268179 
# https://www.medrxiv.org/content/10.1101/2021.12.23.21268179v1 (2021) 
# doi:10.1101/2021.12.23.21268179.

# Last update: 2022-02-16

rm(list=ls())


library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(tidyverse)


  # Import the signature species 
  results.folder <-  "/home/baldanzi/Sleep_apnea/Results/"
  output.plot <- "/home/baldanzi/Sleep_apnea/Results/Plots"
  #wrf <- "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/"

# Function to select the last characters of a string 
  cutlast <- function(char,n){
    l <- nchar(char)
    a <- l-n+1
    return(substr(char,a,l))
  }


  fix.species.name.fun <- function(char){
    char <- gsub("AF22_8LB","AF22-8LB",char)
    char <- gsub("AM42_11","AM42-11",char)
    char <- gsub("TF06_15AC","TF06-15AC",char)
    char <- gsub("AF46_10NS","AF46-10NS",char)
    char <- gsub("AF36_15AT","AF36-15AT",char)
    char <- gsub("_"," ",char)
  }

  # Results from the MGS-AHI/T90/ODI correlation - Full model 

    res.fm <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
    
    Mgs.mgs <- unique(res.fm[,.(MGS,mgs)])

  # Species associated with T90 or ODI 
    
    mgs.t90 <- res.fm[exposure=="t90" & q.value<.05, mgs]
    mgs.odi <- res.fm[exposure=="odi" & q.value<.05, mgs]


  # Import enrichment analysis results GUTSY Atlas 
  ea_gutsy <- fread("/home/baldanzi/Datasets/gutsy_atlas/Supplementary_Table_6.tsv")
  
  names(ea_gutsy) <- gsub("-","_",names(ea_gutsy))
  names(ea_gutsy) <- gsub(" ","_",names(ea_gutsy))
  
  ea_gutsy[,mgs:=paste0("HG3A.",cutlast(Metagenomic_species,4))]
  
  ea_gutsy <- ea_gutsy[mgs %in% unique(c(mgs.t90,mgs.odi)),]
  ea_gutsy <- merge(ea_gutsy, Mgs.mgs, by ="mgs", all.x=T, all.y=F)
  
  
  # Positive and negative associations 
  
  temp.data <- res.fm %>% filter(mgs %in% c(mgs.t90,mgs.odi)) %>% 
    filter(exposure %in% c("t90","odi")) %>%
    select(MGS,mgs,exposure,rho,q.value) %>% 
    pivot_wider(id_cols=c(MGS,mgs), names_from = exposure, values_from = c(rho,q.value))
  
  setDT(temp.data)
  
  # positive species 
  pos.t90 <- temp.data[rho_t90>0 & q.value_t90<.05 & q.value_odi>=.05, MGS]
  pos.t90.odi <- temp.data[rho_t90>0 & q.value_t90<.05 & q.value_odi<.05, MGS]
  pos.odi <- temp.data[rho_odi>0 & q.value_t90>=.05 & q.value_odi<.05, MGS]
  
  pos.mgs <- c(pos.t90,pos.t90.odi,pos.odi)
  
  # negative species 
  neg.t90 <- temp.data[rho_t90<0 & q.value_t90<.05 & q.value_odi>=.05, MGS]
  neg.t90.odi <- temp.data[rho_odi<0 & q.value_t90<.05 & q.value_odi<.05, MGS] # zero
  neg.odi <- temp.data[rho_odi<0 & q.value_t90>=.05 & q.value_odi<.05, MGS]
  
  neg.mgs <- c(neg.t90,neg.t90.odi,neg.odi)
  
  mgs.rel <- c(pos.t90,pos.t90.odi,pos.odi,neg.t90,neg.t90.odi,neg.odi)
  
 
# Enriched pathways/metabolite groups in the POSITIVE correlations ####

  res.table <- ea_gutsy[Direction=="Positive"]
  
  # Filter to pathways/metabolite groups that are FDR associated with at least one signature species 
  
  table.pathways <- res.table %>% group_by(Metabolite_subclass) %>% 
    summarise(nr_p.05 = sum(q_value<.05))
  hm.pathways <- table.pathways$Metabolite_subclass[table.pathways$nr_p.05>=2]

  hm.pathways <- hm.pathways[-which(hm.pathways=="Food Component/Plant")]

  # Matrix results for the heatmap 
  hm.matrix <- res.table %>% select(MGS, Metabolite_subclass, Estimate)%>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = Estimate)
  
  setDF(hm.matrix)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix_pos <- as.matrix(hm.matrix)
  hm.matrix_pos[is.na(hm.matrix_pos)] <- 0
 
  colnames(hm.matrix_pos) <- gsub("__","_", colnames(hm.matrix_pos))
  colnames(hm.matrix_pos) <- gsub("_"," ", colnames(hm.matrix_pos))
  colnames(hm.matrix_pos) <- gsub("Metabolism","Metab.", colnames(hm.matrix_pos))

  # Positive "q-values" to produce the "*" on the heatmap ####
  qvalues_pos <- res.table %>% select(MGS, Metabolite_subclass, q_value)%>% 
      filter(Metabolite_subclass %in% hm.pathways) %>% 
      spread(key=Metabolite_subclass, value = q_value)
  
  setDF(qvalues_pos)

  rownames(qvalues_pos) <- qvalues_pos$MGS
  qvalues_pos$MGS <- NULL

  qvalues_pos <- as.matrix(qvalues_pos)
  qvalues_pos[is.na(qvalues_pos)] <- 1


 # Enriched pathways/metabolite groups in the NEGATIVE correlations ####
 
  res.table <- ea_gutsy[Direction=="Negative"]

  # Filter to pathways/metabolite groups that are FDR associated with at least one signature species 
  table.pathways <- res.table %>% group_by(Metabolite_subclass) %>% 
    summarise(nr_p.05 = sum(q_value<.05))
  hm.pathways <- table.pathways$Metabolite_subclass[table.pathways$nr_p.05>=2]
  
  # Matrix results for the heatmap 
    hm.matrix <- res.table %>% select(MGS, Metabolite_subclass, Estimate)%>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = Estimate)

  setDF(hm.matrix)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix_neg <- as.matrix(hm.matrix)
  hm.matrix_neg[is.na(hm.matrix_neg)] <- 0

  colnames(hm.matrix_neg) <- gsub("__","_",colnames(hm.matrix_neg))
  colnames(hm.matrix_neg) <- gsub("_"," ",colnames(hm.matrix_neg))
  colnames(hm.matrix_neg) <- gsub("Metabolism","Metab.",colnames(hm.matrix_neg))
  colnames(hm.matrix_neg) <- gsub(" (PC)","",colnames(hm.matrix_neg))

  # Negative "q-values" to produce the "*" on the heatmap ####
    qvalues_neg <- res.table %>% select(MGS, Metabolite_subclass, q_value) %>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = q_value)
  
  setDF(qvalues_neg)

  rownames(qvalues_neg) <- qvalues_neg$MGS
  qvalues_neg$MGS <- NULL

  qvalues_neg <- as.matrix(qvalues_neg)
    qvalues_neg[is.na(qvalues_neg)] <- 1


# Row Annotation ####

  mgs.rel <- mgs.rel[mgs.rel %in% rownames(hm.matrix)]

  annotation <- data.table(mgs.rel = mgs.rel, 
                         correlation = vector(mode="character", 
                                             length=length(mgs.rel)))
  
  list.mgs <- list(pos.t90=pos.t90,pos.t90.odi=pos.t90.odi,pos.odi=pos.odi,
                   neg.t90=neg.t90,neg.t90.odi=neg.t90.odi,neg.odi=neg.odi)
 
  for(i in 1:length(list.mgs)){
    annotation[mgs.rel %in% list.mgs[[i]], correlation := names(list.mgs)[i]]
    }
  
  annotation[,correlation:=factor(correlation)]
  
  setDF(annotation)
  
  rownames(annotation) <- annotation$mgs.rel 

  # Assert that the heatmap matrix and the annotation vector have in the same order
  hm.matrix_neg <- hm.matrix_neg[match(annotation$mgs.rel,rownames(hm.matrix_neg)), ]
  hm.matrix_pos <- hm.matrix_pos[match(annotation$mgs.rel,rownames(hm.matrix_pos)), ]

  # Assert that the q-values matrix and the annotation vector have in the same order
  qvalues_pos <- qvalues_pos[match(annotation$mgs.rel,rownames(qvalues_pos)),]
  qvalues_neg <- qvalues_neg[match(annotation$mgs.rel,rownames(qvalues_neg)),]

  map.row = c("cornflowerblue","violetred4","grey83")
  names(map.row) = substr(names(list.mgs),5,20)[1:3]

  # Fix species names to better visualization in the heatmap
  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(hm.matrix_pos)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(hm.matrix_pos) <- noms

  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(hm.matrix_neg)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(hm.matrix_neg) <- noms

  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(annotation)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(annotation) <- noms

  
  
  # Create object with the annotations for the Heatmap 
  value=factor(substr(annotation$correlation,5,20),levels=c("t90","t90.odi","odi"))
  
  ha = HeatmapAnnotation(#empty=anno_empty(border=F, width = unit(6,"mm")),
                         value=value, which="row",
                         simple_anno_size = unit(2,"mm"),
                         col =  list(value=map.row),
                         show_annotation_name = F,
                         annotation_legend_param = list(title="",
                                                        labels = c("T90","Both","ODI")))
  
      # Create heatmap
  set.seed(2)
  ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
  h1 = Heatmap(hm.matrix_pos, 
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(qvalues_pos[i, j] < 0.05) {
                 grid.text('*', x, y)
               } },
             column_names_rot =  40, 
             column_names_side = "bottom",
             cluster_columns = T,
             show_column_dend = F,
             
             #row_split = factor(annotation$correlation),
             
             row_split = factor(c(rep("Positive T90/ODI-associated species",length(pos.mgs)),
                                     rep("Negative T90/ODI-associated species",length(neg.mgs)))),
             cluster_row_slices = T,
             cluster_rows = T,
             #row_title = c("Positive T90/ODI-associated species",
                         #  "Negative T90/ODI-associated species"),
             row_title_gp = gpar(fontsize=10, fill="white",border="black"),
             row_gap = unit(3, "mm"),
             
              row_names_side = "left",
             show_row_dend = F,
             
             col = circlize::colorRamp2(c(0,1.5,3),c("white","indianred2", "red3")),
             name="NES (pos)",
             column_labels = colnames(hm.matrix_pos),
             column_names_gp = gpar(fontsize = 8),
             row_names_gp = gpar(fontsize = 8),
             column_title = "Positive\nspecies-metabolites\nassociations",
             column_title_gp = gpar(fontsize = 9,fontface='bold'),
             heatmap_legend_param = list(labels_gp = gpar(fontsize=7),
                                         title_gp = gpar(fontsize=8),
                                         grid_width = unit(1.5, "mm"))  ,
             left_annotation = ha, 
            width = unit(5*ncol(hm.matrix_pos),"mm")) + # The symbol '%v%' means one heatmap over the other
  
  Heatmap(hm.matrix_neg, 
          cell_fun = function(j, i, x, y, w, h, fill) {
            if(qvalues_neg[i, j] < 0.05) {
              grid.text('*', x, y)
            } },
          column_names_rot =  40,
          column_names_side = "bottom",
          cluster_columns = T,
          show_column_dend = F,
          
          row_names_side = "left",
          show_row_dend = F,
          
          col = circlize::colorRamp2(c(0,1.5,3),c("white","dodgerblue1","blue4" )),
          name="NES (neg)",
          column_labels = colnames(hm.matrix_neg),
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 8),
          
          column_title = "Negative\nspecies-metabolites\nassociations",
          column_title_gp = gpar(fontsize = 9,fontface='bold'),
          heatmap_legend_param = list(labels_gp = gpar(fontsize=7),
                                      title_gp = gpar(fontsize=8),
                                      grid_width = unit(1.5, "mm")),
          width = unit(5*ncol(hm.matrix_neg),"mm")) 

 

  # Save the plot in pdf 
  pdf(file = paste0(output.plot, "mgs_subpathway_ea_heatmap_gutsy_new.pdf"), 
    width = 10, height = 10)
  set.seed(2)
  draw(h1, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend=T)
  dev.off()

  # Save the plot in png 
  png(filename =  "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/mgs_subpathway_ea_heatmap_gutsy.png", 
  width = 9, height = 8, units = 'in', res = 1500 )
  set.seed(1)
  draw(h1, 
     column_title_gp = gpar(fontsize = 12, fontface="bold"),
     heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend=T)

  dev.off()


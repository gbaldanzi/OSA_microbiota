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
  wrf <- "/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512/"

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
    char <- gsub("Marseille_P4005","Marseille-P4005", char)
    char <- gsub('[(]HG3A.0168[)]',"",char)
    char <- gsub("_"," ",char)
  }

  # Results from the MGS-AHI/T90/ODI correlation - Full model 

    res.bm <- fread(paste0(results.folder,"cor.bmi_all.var_mgs.tsv"))
    
    Mgs.mgs <- unique(res.bm[,.(MGS,mgs)])

  # Species associated with T90 or ODI 
    
    mgs.t90 <- res.bm[exposure=="t90" & q.value<.05, mgs]
    mgs.odi <- res.bm[exposure=="odi" & q.value<.05, mgs]


  # Import enrichment analysis results GUTSY Atlas 
  ea_gutsy <- fread("/home/baldanzi/Datasets/gutsy_atlas/Supplementary_Table_6.tsv")
  
  names(ea_gutsy) <- gsub("-","_",names(ea_gutsy))
  names(ea_gutsy) <- gsub(" ","_",names(ea_gutsy))
  
  ea_gutsy[,mgs:=paste0("HG3A.",cutlast(Metagenomic_species,4))]
  
  ea_gutsy <- ea_gutsy[mgs %in% unique(c(mgs.t90,mgs.odi)),]
  ea_gutsy <- merge(ea_gutsy, Mgs.mgs, by ="mgs", all.x=T, all.y=F)
  
  
  # Positive and negative associations 
  
  temp.data <- res.bm[,.(MGS,mgs,exposure,rho,q.value)] %>% 
    filter(mgs %in% c(mgs.t90,mgs.odi)) %>% 
    filter(exposure %in% c("t90","odi")) %>%
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
  hm.pathways <- table.pathways$Metabolite_subclass[table.pathways$nr_p.05>=1]

  #hm.pathways <- hm.pathways[-which(hm.pathways=="Food Component/Plant")]

  # Matrix results for the heatmap 
  hm.matrix <- res.table[,.(MGS, Metabolite_subclass, Estimate)] %>% 
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
  qvalues_pos <- res.table[,.(MGS, Metabolite_subclass, q_value)] %>% 
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
  hm.pathways <- table.pathways$Metabolite_subclass[table.pathways$nr_p.05>=1]
  
  # Matrix results for the heatmap 
    hm.matrix <- res.table[,.(MGS, Metabolite_subclass, Estimate)] %>% 
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
    qvalues_neg <- res.table[,.(MGS, Metabolite_subclass, q_value)] %>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = q_value)
  
  setDF(qvalues_neg)

  rownames(qvalues_neg) <- qvalues_neg$MGS
  qvalues_neg$MGS <- NULL

  qvalues_neg <- as.matrix(qvalues_neg)
    qvalues_neg[is.na(qvalues_neg)] <- 1

    
    # Association with T90 and ODI 
    
    res.bm[,exposure:= factor(exposure, levels = c("t90","odi","ahi"))]
    
     t90.odi.cor <- res.bm[MGS %in% mgs.rel & exposure %in% c("t90","odi"), 
                           .(MGS,exposure,rho)]
     t90.odi.q <- res.bm[MGS %in% mgs.rel & exposure %in% c("t90","odi"), 
                           .(MGS,exposure,q.value)]
     
     t90.odi.cor[,exposure := factor(toupper(exposure),levels=c("T90","ODI"))]
     t90.odi.q[,exposure := factor(toupper(exposure),levels=c("T90","ODI"))]
     
     t90.odi.cor <- t90.odi.cor %>% spread(key=exposure, value = rho)
     t90.odi.q <- t90.odi.q %>% spread(key=exposure, value = q.value)
     
     setDT(t90.odi.cor)
     setDT(t90.odi.q)
     
     mgs.t90.odi.cor <- t90.odi.cor$MGS 
     mgs.t90.odi.q <- t90.odi.q$MGS
     
     t90.odi.cor$MGS <- t90.odi.q$MGS <-  NULL
     
     t90.odi.cor <- as.matrix(t90.odi.cor)
     t90.odi.q <- as.matrix(t90.odi.q)
     
    rownames(t90.odi.cor) <- mgs.t90.odi.cor
    rownames(t90.odi.q) <- mgs.t90.odi.q

  # Assert that the heatmap matrixes have the same order 
  hm.matrix_neg <- hm.matrix_neg[match(mgs.rel,rownames(hm.matrix_neg)), ]
  hm.matrix_pos <- hm.matrix_pos[match(mgs.rel,rownames(hm.matrix_pos)), ]

  qvalues_pos <- qvalues_pos[match(mgs.rel,rownames(qvalues_pos)),]
  qvalues_neg <- qvalues_neg[match(mgs.rel,rownames(qvalues_neg)),]
  
  t90.odi.cor <- t90.odi.cor[match(mgs.rel, rownames(t90.odi.cor)),]
  t90.odi.q <- t90.odi.q[match(mgs.rel, rownames(t90.odi.q)),]


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

  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(t90.odi.cor)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(t90.odi.cor) <- noms

  
  
      # Create heatmap
  set.seed(10)
  ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
  h1 = Heatmap(t90.odi.cor, 
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(t90.odi.q[i, j] < 0.05) {
                   grid.text('*', x, y)
                 } },
               column_names_rot=40, 
               column_names_side = "bottom",
               cluster_columns = FALSE, 
               show_column_dend = F, 
               column_names_gp = gpar(fontsize = 9),
               
               row_split = factor(c(rep("Positive T90/ODI-associated species",length(pos.mgs)),
                                    rep("Negative T90/ODI-associated species",length(neg.mgs))), 
                                  levels = c("Positive T90/ODI-associated species", 
                                             "Negative T90/ODI-associated species")),
               cluster_row_slices = F,
               cluster_rows = T,   
               row_title_gp = gpar(fontsize=10, fill="white",border="black"),
               row_gap = unit(3, "mm"),
               
               col = circlize::colorRamp2( c(-0.1, 0, 0.1), c("#482173FF","white","goldenrod2")),
               #name = expression(rho),
               
               heatmap_legend_param = list( title = expression(rho),
                                          labels_gp = gpar(fontsize=7),
                                           title_gp = gpar(fontsize=10),
                                           grid_width = unit(3, "mm"), 
                                           direction = "vertical"),
               
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 8),
               show_row_dend = F, 
               width = unit(5*ncol(t90.odi.cor),"mm")) +
    
        Heatmap(hm.matrix_pos, 
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(qvalues_pos[i, j] < 0.05) {
                 grid.text('*', x, y)
               } },
             column_names_rot =  40, 
             column_names_side = "bottom",
             cluster_columns = T,
             show_column_dend = F,

             #row_title = c("Positive T90/ODI-associated species",
                         #  "Negative T90/ODI-associated species"),
             row_title_gp = gpar(fontsize=10, fill="white",border="black"),
             row_gap = unit(3, "mm"),
             
             row_names_side = "left",

             col = circlize::colorRamp2(c(3,1.5,0),c("red3","indianred2","white")),
             name="NES (pos)",
             column_labels = colnames(hm.matrix_pos),
             column_names_gp = gpar(fontsize = 8),
             row_names_gp = gpar(fontsize = 8),
             column_title = "Positive",
             column_title_gp = gpar(fontsize = 10,fontface='bold'),
             heatmap_legend_param = list( title = "NES",
                                         labels_gp = gpar(fontsize=7),
                                         title_gp = gpar(fontsize=9, fontface="bold"),
                                         grid_width = unit(2, "mm"), 
                                         direction = "horizontal")  ,
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
          
          col = circlize::colorRamp2(c(0,1.5,3),c("white","dodgerblue1","blue4" )),
          name="NES (neg)",
          column_labels = colnames(hm.matrix_neg),
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 8),
          row_names_side = "left",
          
          column_title = "Negative",
          column_title_gp = gpar(fontsize = 10,fontface='bold'),
          heatmap_legend_param = list(title = "",
                                      labels_gp = gpar(fontsize=7),
                                      title_gp = gpar(fontsize=8),
                                      grid_width = unit(2, "mm"), 
                                      direction = "horizontal"),
          width = unit(5*ncol(hm.matrix_neg),"mm")) 

  

  message("Saving plots")
  # Save the plot in pdf 
  pdf(file = paste0(output.plot, "mgs_subpathway_ea_heatmap_gutsy.pdf"), 
    width = 10, height = 10)
  set.seed(10)
  draw(h1, heatmap_legend_side = "right", merge_legend=F)
  dev.off()
  
  pdf(file = paste0(wrf, "mgs_subpathway_ea_heatmap_gutsy.pdf"), 
      width = 11, height = 9)
  set.seed(10)
  draw(h1, heatmap_legend_side = "right", merge_legend=F, column_title = "Species-metabolites associations:")
  dev.off()

  # Save the plot in png 
  #png(filename =  "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/mgs_subpathway_ea_heatmap_gutsy.png", 
  #width = 9, height = 8, units = 'in', res = 1500 )
  #set.seed(1)
  #draw(h1, 
  #   column_title_gp = gpar(fontsize = 12, fontface="bold"),
  #   heatmap_legend_side = "top", merge_legend=F)

  #dev.off()
  message("Plots saved")


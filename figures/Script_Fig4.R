# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script will create the Heatmap (Fig 4) containing the enrichment
# analysis results 


# The species enrichment analysis results are imported from the Gutsy Atlas 
# reference: Dekkers, K. F. et al. An online atlas of human plasma metabolite 
# signatures of gut microbiome composition. 2021.12.23.21268179 
# https://www.medrxiv.org/content/10.1101/2021.12.23.21268179v1 (2021) 
# doi:10.1101/2021.12.23.21268179.


rm(list=ls())


library(RColorBrewer)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(data.table)
library(tidyverse)
library(scales)
library(cowplot)


  # Import the signature species 
  #results.folder <- "/Users/gabba126/Documents/PhD_projects/1.Sleep/Results/"
  results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
  output.plot <- "/home/baldanzi/Sleep_apnea/Results/Plots/"
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
  
    fix.subclass.name.fun <- function(char){
      char <- gsub("__PC_", "",char)
      char <- gsub("___","_",char)
      char <- gsub("__","_",char)
      char <- gsub("Drug", "Drug -", char)
      char <- gsub("_"," ",char)
      char <- gsub(" PI ", " (PI)",char)
      char <- gsub(" PE ", " (PE)",char)
      char <- gsub("Analgesics","Analgesics,", char)
      return(char)
    }
      
  # Taxonomy 
    taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy")

  # Results from the MGS-AHI/T90/ODI correlation - extended model 

    res.bm <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
    
    #res.bm <- merge(res.bm,taxonomy,by.x="MGS", by.y = "maintax_mgs", all.x = T, all.y = F)
    
    res.bm[,mgs := cutlast(MGS,9)]
    
    Mgs.mgs <- unique(res.bm[,.(MGS,mgs)])

  # Species in association with T90 or ODI (any direction)
    sig.mgs <- unique(res.bm[exposure %in% c("t90","odi") & q.value<.05, mgs])
    
    
  # Species in negative association with T90 or ODI 
    neg.mgs <- unique(res.bm[exposure %in% c("t90","odi") & rho<0 & q.value<.05, mgs])
    
    
  # Species in negative association with T90 or ODI 
    pos.mgs <- unique(res.bm[exposure %in% c("t90","odi") & rho>0 & q.value<.05, mgs])
  

  # Import enrichment analysis results GUTSY Atlas 
  ea_gutsy <- fread(paste0(results.folder,"Supplementary_Table_6.tsv"))
  
  names(ea_gutsy) <- gsub("-","_",names(ea_gutsy))
  names(ea_gutsy) <- gsub(" ","_",names(ea_gutsy))
  
  ea_gutsy[,mgs:=paste0("HG3A.",cutlast(Metagenomic_species,4))]
  
  # Restrict the GUTSY Atlas table to the species of interest 
  ea_gutsy <- ea_gutsy[mgs %in% sig.mgs,]
  ea_gutsy <- merge(ea_gutsy, Mgs.mgs, by ="mgs", all.x=T, all.y=F)
  
  # Species in NEGATIVE association with T90/ODI ####

  # Enriched pathways/metabolite groups in the POSITIVE correlations ####
  res.table <- ea_gutsy[Direction=="Positive" & mgs %in% neg.mgs, .(MGS, Metabolite_subclass, Estimate, q_value)]
  
  # Filter to pathways/metabolite groups that are FDR associated with at least two  species 
  
  res.table.sig <- res.table[q_value<.05, ] %>% group_by(Metabolite_subclass) %>% count() %>% filter(n>=2) %>% select(Metabolite_subclass)
  
  hm.pathways <- res.table.sig$Metabolite_subclass
  hm.pathways <- hm.pathways[-which(hm.pathways == "Food Component/Plant")]


  # Matrix results for the heatmap 
  hm.matrix <- res.table[,.(MGS, Metabolite_subclass, Estimate)] %>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = Estimate)
  
  setDF(hm.matrix)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix_pos2 <- as.matrix(hm.matrix)
  hm.matrix_pos2[is.na(hm.matrix_pos2)] <- 0
 
  colnames(hm.matrix_pos2) <- gsub("__","_", colnames(hm.matrix_pos2))
  colnames(hm.matrix_pos2) <- gsub("_"," ", colnames(hm.matrix_pos2))
  colnames(hm.matrix_pos2) <- gsub("Metabolism","Metab.", colnames(hm.matrix_pos2))

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
 
  res.table <- ea_gutsy[Direction=="Negative" & mgs %in% neg.mgs, .(MGS, Metabolite_subclass, Estimate, q_value)]
  
  # Filter to pathways/metabolite groups that are FDR associated with at least two  species
  res.table.sig <- res.table[q_value<.05, ] %>% group_by(Metabolite_subclass) %>% count() %>% filter(n>2) %>% select(Metabolite_subclass)
  
  hm.pathways <- res.table.sig$Metabolite_subclass
  
  # Matrix results for the heatmap 
    hm.matrix <- res.table[,.(MGS, Metabolite_subclass, Estimate)] %>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = Estimate)

  setDF(hm.matrix)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix_neg2 <- as.matrix(hm.matrix)
  hm.matrix_neg2[is.na(hm.matrix_neg2)] <- 0

  colnames(hm.matrix_neg2) <- gsub("__","_",colnames(hm.matrix_neg2))
  colnames(hm.matrix_neg2) <- gsub("_"," ",colnames(hm.matrix_neg2))
  colnames(hm.matrix_neg2) <- gsub("Metabolism","Metab.",colnames(hm.matrix_neg2))
  colnames(hm.matrix_neg2) <- gsub(" (PC)","",colnames(hm.matrix_neg2))

  # Negative "q-values" to produce the "*" on the heatmap ####
    qvalues_neg <- res.table[,.(MGS, Metabolite_subclass, q_value)] %>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = q_value)
  
    setDF(qvalues_neg)

    rownames(qvalues_neg) <- qvalues_neg$MGS
    qvalues_neg$MGS <- NULL

    qvalues_neg <- as.matrix(qvalues_neg)
    qvalues_neg[is.na(qvalues_neg)] <- 1

    
    
  # Assert that the heatmap matrixes have the same row order 
    
    names.neg.mgs <-  Mgs.mgs[mgs %in% neg.mgs, MGS]
    
    assert.order <- function(x){
          return(x[match(names.neg.mgs, rownames(x)),])
          }
    
    
    hm.matrix_pos2 <- assert.order(hm.matrix_pos2)
    hm.matrix_neg2 <- assert.order(hm.matrix_neg2) 
  
    qvalues_pos <- assert.order(qvalues_pos)
    qvalues_neg <- assert.order(qvalues_neg)
  
  
  # Merge hm.matrix_pos and hm.matrix_neg
  merged.matrix <- cbind(hm.matrix_pos2,(hm.matrix_neg2*-1))
  merged.matrix.q <- cbind(qvalues_pos,qvalues_neg)
  
  
  # Fix species names to better visualization in the heatmap
  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(merged.matrix)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(merged.matrix) <- noms
  
  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(merged.matrix.q)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(merged.matrix.q) <- noms
  
  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , taxonomy$maintax_mgs)))
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  taxonomy$name <- noms
  

  # Only include species with at least one FDR sig association
  species.sig <- merged.matrix.q<0.05
  species.sig <- apply(species.sig, 1, any)==T
  species.sig <- names(species.sig)[species.sig==T]
  
  merged.matrix2 <- merged.matrix[rownames(merged.matrix) %in% species.sig,]
  merged.matrix.q2 <- merged.matrix.q[rownames(merged.matrix.q) %in% species.sig,]
  
  
  # Annotation #### 
  
  taxa <- taxonomy[name %in% rownames(merged.matrix2), .(name,family,order,class,phylum)]
  order2 <- taxa[match(rownames(merged.matrix2),name ), order]
  ha2 <- HeatmapAnnotation(which = "column",  Order = order2, simple_anno_size = unit(2,"mm"),
                           col= list(Order = c("Bacteroidales" = "#1B9E77" ,
                                     "Candidatus Borkfalkiales" = "#D95F02", 
                                     "Eubacteriales" = "#66A61E",   
                                     "Eggerthellales" = "#E7298A", 
                                     "Erysipelotrichales" = "#7570B3",
                                     "Victivallales" = "#E6AB02", 
                                     "unclassified" = "gray80")), 
                           show_legend = F, 
                           annotation_name_side = "left", 
                           annotation_name_gp = gpar(fontsize=8))
  
  

  # Create heatmap ####
  set.seed(10)
  
  ht_opt$TITLE_PADDING = unit(c(2, 2), "mm")
  ht_opt$COLUMN_ANNO_PADDING = unit(1, "mm")
  ht_opt$DIMNAME_PADDING = unit(1, "mm")
  
  h2 <-   Heatmap(t(merged.matrix2), 
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(t(merged.matrix.q2)[i, j] < 0.05) {
                 grid.text('*', x, y)
               } },
             # rownq =  0, 
             # row = "bottom",
             cluster_columns = T,
             show_column_dend = F,
             show_row_names = T,
             show_column_names = F,
             show_row_dend = FALSE,
             cluster_rows = TRUE,
             cluster_row_slices = F,
             top_annotation = ha2,
             
             row_split = factor(c(rep("pos", ncol(hm.matrix_pos2)),
                                  rep("neg",ncol(hm.matrix_neg2))), 
                                levels = c("pos","neg")),
             row_gap = unit(2, "mm"),
             #right_annotation = ha1n,

             #row_title = c("positive\nassociations", "negative\nassociations"),
                         #  "Negative T90/ODI-associated species"),
            
             #row_gap = unit(3, "mm"),
             
             row_names_side = "left",

             col = circlize::colorRamp2(c(-2.5,-1.6,0,1.6,2.5),c(muted("blue4"), "dodgerblue1", "white",
                                                             "indianred2",muted("red3"))),
             # name="NES (pos)",
             column_labels = colnames(merged.matrix2),
             column_names_gp = gpar(fontsize = 8),
             row_names_gp = gpar(fontsize = 8),
             column_title = "Species in negative association with T90/ODI",
             row_title = c("positive","negative"),
             row_title_gp = gpar(fontsize=8, border = "black"),
             column_title_gp = gpar(fontsize = 9, fontface="bold" ),
             
             #heatmap_legend_param = list( title = "NES",
                                        # labels_gp = gpar(fontsize=7),
                                        # title_gp = gpar(fontsize=9, fontface="bold"),
                                        # grid_width = unit(2, "mm"), 
                                        # direction = "horizontal")  ,
             width = unit(1.4*nrow(merged.matrix2),"mm"), 
             height = unit(3.5*ncol(merged.matrix2),"mm"))  
  

  
  
  
  
  # Species in POSITIVE association with T90/ODI ####
  

  # Enriched pathways/metabolite groups in the POSITIVE correlations ####
  res.table <- ea_gutsy[Direction=="Positive" & mgs %in% pos.mgs, .(MGS, Metabolite_subclass, Estimate, q_value)]
  

  # Filter to pathways/metabolite groups that are FDR associated with at least one signature species 
  res.table.sig <- res.table[q_value<.05, ] %>% group_by(Metabolite_subclass) %>% count() %>% filter(n>=1) %>% select(Metabolite_subclass)
  
  hm.pathways <- res.table.sig$Metabolite_subclass
  
  hm.pathways <- hm.pathways[-(grep("Drug ",hm.pathways))]
  
  
  # Matrix results for the heatmap 
  hm.matrix <- res.table[,.(MGS, Metabolite_subclass, Estimate)] %>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = Estimate)
  
  setDF(hm.matrix)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix_pos1 <- as.matrix(hm.matrix)
  hm.matrix_pos1[is.na(hm.matrix_pos1)] <- 0
  
  colnames(hm.matrix_pos1) <- gsub("__","_", colnames(hm.matrix_pos1))
  colnames(hm.matrix_pos1) <- gsub("_"," ", colnames(hm.matrix_pos1))
  colnames(hm.matrix_pos1) <- gsub("Metabolism","Metab.", colnames(hm.matrix_pos1))
  
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
  
  res.table <- ea_gutsy[Direction=="Negative" & mgs %in% pos.mgs, .(MGS, Metabolite_subclass, Estimate, q_value)]
  

  # Filter to pathways/metabolite groups that are FDR associated with at least one signature species 
  hm.pathways <- unique(res.table[q_value<.05, Metabolite_subclass])
  
  hm.pathways <- hm.pathways[-(grep("Drug ",hm.pathways))]
  
  # Matrix results for the heatmap 
  hm.matrix <- res.table[,.(MGS, Metabolite_subclass, Estimate)] %>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = Estimate)
  
  setDF(hm.matrix)
  
  rownames(hm.matrix) <- hm.matrix$MGS
  hm.matrix$MGS <- NULL
  
  hm.matrix_neg1 <- as.matrix(hm.matrix)
  hm.matrix_neg1[is.na(hm.matrix_neg1)] <- 0
  
  colnames(hm.matrix_neg1) <- gsub("__","_",colnames(hm.matrix_neg1))
  colnames(hm.matrix_neg1) <- gsub("_"," ",colnames(hm.matrix_neg1))
  colnames(hm.matrix_neg1) <- gsub("Metabolism","Metab.",colnames(hm.matrix_neg1))
  colnames(hm.matrix_neg1) <- gsub(" (PC)","",colnames(hm.matrix_neg1))
  
  # Negative "q-values" to produce the "*" on the heatmap ####
  qvalues_neg <- res.table[,.(MGS, Metabolite_subclass, q_value)] %>% 
    filter(Metabolite_subclass %in% hm.pathways) %>% 
    spread(key=Metabolite_subclass, value = q_value)
  
  setDF(qvalues_neg)
  
  rownames(qvalues_neg) <- qvalues_neg$MGS
  qvalues_neg$MGS <- NULL
  
  qvalues_neg <- as.matrix(qvalues_neg)
  qvalues_neg[is.na(qvalues_neg)] <- 1
  
  

  
  # Assert that the heatmap matrixes have the same order 
  names.pos.mgs <-  Mgs.mgs[mgs %in% pos.mgs, MGS]
  
  assert.order <- function(x){
    return(x[match(names.pos.mgs, rownames(x)),])
  }
  
  hm.matrix_pos1 <- assert.order(hm.matrix_pos1)
  hm.matrix_neg1 <- assert.order(hm.matrix_neg1) 
  
  qvalues_pos <- assert.order(qvalues_pos)
  qvalues_neg <- assert.order(qvalues_neg)
  
  
  # Merge hm.matrix_pos1 and hm.matrix_neg1
  merged.matrix1 <- cbind(hm.matrix_pos1,(hm.matrix_neg1*-1))
  merged.matrix.q1 <- cbind(qvalues_pos,qvalues_neg)
  
  
  
  # Fix species names to better visualization in the heatmap
  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(merged.matrix1)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(merged.matrix1) <- noms
  
  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(merged.matrix.q1)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(merged.matrix.q1) <- noms
  
  
  # Annotation ####
  
  taxa <- taxonomy[name %in% rownames(merged.matrix1), .(name,family,order,class,phylum)]
  order1 <- taxa[match(rownames(merged.matrix1),name ), order]
  ha1 <- HeatmapAnnotation(which = "column", Order = order1, simple_anno_size = unit(2,"mm"),
            col= list( Order = c("Coriobacteriales" = "#A6761D" ,
            "Bacillales" = "#666666", 
            "Eubacteriales" = "#66A61E",
            "Erysipelotrichales" = "#7570B3",  
            "Lactobacillales" = "#386CB0")) , 
            show_legend = F, 
            annotation_name_side = "left", 
            annotation_name_gp = gpar(fontsize = 8))
  

 
  # Create heatmap ####
  set.seed(10)
  
  ht_opt$TITLE_PADDING = unit(c(2, 2), "mm")
  ht_opt$COLUMN_ANNO_PADDING = unit(1, "mm")
  ht_opt$DIMNAME_PADDING = unit(1, "mm")
  
  h1 <-   Heatmap(t(merged.matrix1), 
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(t(merged.matrix.q1)[i, j] < 0.05) {
                    grid.text('*', x, y)
                  } },
                # rownq =  0, 
                # row = "bottom",
                cluster_columns = T,
                show_column_dend = F,
                show_row_names = T,
                show_column_names = F,
                show_row_dend = FALSE,
                cluster_rows = TRUE,
                cluster_row_slices = F,
                
                row_split = factor(c(rep("pos", ncol(hm.matrix_pos1)),
                                     rep("neg",ncol(hm.matrix_neg1))), 
                                   levels = c("pos","neg")),
                row_gap = unit(2, "mm"),

                row_names_side = "left",
                
                top_annotation = ha1, 
                
                col = circlize::colorRamp2(c(-2.5,-1.6,0,1.6,2.5),c(muted("blue4"), "dodgerblue1", "white",
                                                                    "indianred2",muted("red3"))),
                # name="NES (pos)",
                column_labels = colnames(merged.matrix2),
                column_names_gp = gpar(fontsize = 8),
                row_names_gp = gpar(fontsize = 8),
                column_title = "Species in positive\nassociation with T90/ODI",
                row_title = c("positive","negative"),
                row_title_gp = gpar(fontsize=8, border = "black"),
                column_title_gp = gpar(fontsize = 9, fontface="bold" ),
                
                #heatmap_legend_param = list( title = "NES",
                # labels_gp = gpar(fontsize=7),
                # title_gp = gpar(fontsize=9, fontface="bold"),
                # grid_width = unit(2, "mm"), 
                # direction = "horizontal")  ,
                width = unit(1.4*nrow(merged.matrix1),"mm"), 
                height = unit(3.5*ncol(merged.matrix1),"mm"))  
  
  
  
  
  # Legend ####
  col_fun = circlize::colorRamp2(c(0,1.6,2.5),c("white","indianred2", muted("red3")))
  col_fun2 = circlize::colorRamp2(c(0,1.6,2.5),c("white","dodgerblue1", muted("blue4")))
  #col_fun3 = circlize::colorRamp2(c(-0.1, 0, 0.1), c(muted("darkcyan"),"white","goldenrod2"))
  
  
  lgd = Legend(col_fun = col_fun, title = "NES", at = c(0,1,2,3), direction = "horizontal",
               labels_gp = gpar(fontsize=6), title_gp = gpar(fontsize=8), 
               grid_height = unit(2.5,"mm"))
  lgd2 = Legend(col_fun = col_fun2, at = c(0,1,2,3), direction = "horizontal",
                labels_gp = gpar(fontsize=6), 
                grid_height = unit(2.5,"mm"))
  
  orders <- unique(c(order1,order2))
  orders <- orders[order(orders)]
  
  lgd3 = Legend(labels = orders, title = "Order", labels_gp = gpar(fontsize = 6.5), title_gp = gpar(fontsize=8),
                legend_gp = gpar(fill = c( '#666666', '#1B9E77' , "#D95F02", "#A6761D",  "#E7298A" , 
                                            "#7570B3" ,"#66A61E", '#386CB0' , "gray80" , "#E6AB02")),
                grid_height = unit(2.5,"mm") , grid_width = unit(2.5,"mm"))
 
  
  #lgd3 = Legend(col_fun = col_fun3, at = c(-.1,-0.05,0,0.05,0.1), direction = "horizontal", 
  #           title = expression(rho), legend_height = unit(20,"mm"), labels_gp = gpar(fontsize=7))
  
  pd1 = packLegend(lgd, lgd2, direction = "vertical")
  pd2 = packLegend(pd1, lgd3, direction = "horizontal", column_gap = unit(7, "mm"))
  
  w = convertWidth(unit(1, "npc")*(9/10), "mm", valueOnly = TRUE)
  h = convertHeight(unit(1, "npc")*(4/5), "mm", valueOnly = TRUE)
  
  p1 = grid.grabExpr(draw(h1, show_heatmap_legend = F, row_title = "Species-metabolite associations         ", 
                          row_title_gp = gpar(fontsize=8), 
                          padding = unit(c(3,.1,.1,32),"mm")), 
                     width = w, height = h)
  
  p2 = grid.grabExpr(draw(h2, show_heatmap_legend=F, row_title = "Species-metabolite associations           ", 
                     row_title_gp = gpar(fontsize=8), padding = unit(c(2,.15,.1,.1),"mm")),
                     width = unit(10,"mm"), height = unit(6,"mm"))
  
  p3 = grid.grabExpr(draw(pd2, x = unit(-.35, "npc"), y = unit(.5, "npc"), just = c("left")))
                     #padding = unit(c(3,.1,.1,5),"mm")))
  
  pfinal <- plot_grid(p2,plot_grid(p1,p3,nrow=1, rel_widths = c(1,.2), labels = c("b","")),
                      ncol=1, labels = c("a",""), label_y = .94)
  
  
  save_plot("Fig4.pdf", pfinal, base_height = 4.67)
 

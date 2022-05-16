# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script will create the Heatmap (Fig 3) containing the enrichment
# analysis results and the associations with health outcomes 

# On the "x-axis", we display the metabolite groups and systolic blood pressure,
# diastolic blood pressure and Hb1Ac

# On the "y-axis", we display the species and the microbial pathways. 

# The species enrichment analysis results are imported from the Gutsy Atlas 
# reference: Dekkers, K. F. et al. An online atlas of human plasma metabolite 
# signatures of gut microbiome composition. 2021.12.23.21268179 
# https://www.medrxiv.org/content/10.1101/2021.12.23.21268179v1 (2021) 
# doi:10.1101/2021.12.23.21268179.


rm(list=ls())


library(RColorBrewer)
library(ComplexHeatmap)
suppressPackageStartupMessages(library(circlize))
library(data.table)
library(tidyverse)
library(scales)


  # Import the signature species 
  results.folder <- "/Users/gabba126/Documents/PhD_projects/1.Sleep/Results/"
  # results.folder <-  '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/'
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
      

  # Results from the MGS-AHI/T90/ODI correlation - extended model 

    res.bm <- fread(paste0(results.folder,"cor2_all.var_mgs.tsv"))
    
    res.bm[,mgs := cutlast(MGS,9)]
    
    Mgs.mgs <- unique(res.bm[,.(MGS,mgs)])

  # Species associated with T90 or ODI 
    
    mgs.t90 <- res.bm[exposure=="t90" & rho>0 & q.value<.05, mgs]
    mgs.odi <- res.bm[exposure=="odi" & rho>0 & q.value<.05, mgs]


  # Import enrichment analysis results GUTSY Atlas 
  ea_gutsy <- fread(paste0(results.folder,"Supplementary_Table_6.tsv"))
  
  names(ea_gutsy) <- gsub("-","_",names(ea_gutsy))
  names(ea_gutsy) <- gsub(" ","_",names(ea_gutsy))
  
  ea_gutsy[,mgs:=paste0("HG3A.",cutlast(Metagenomic_species,4))]
  
  # Restrict the GUTSY Atlas table to the species of interest 
  ea_gutsy <- ea_gutsy[mgs %in% unique(c(mgs.t90,mgs.odi)),]
  ea_gutsy <- merge(ea_gutsy, Mgs.mgs, by ="mgs", all.x=T, all.y=F)
  
  # Import results on the enrichment of metabolites subclasses ####
  # ea.gmm.met <- fread(paste0(results.folder, "ea_gmm_subclass.tsv"))
  
  # gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  # setnames(gmm.names,"Module","modules")
  # gmm.names[,Name:=str_to_title(Name)]
  # gmm.names[,Name:=gsub("Ii","II",Name)]
  
  # ea.gmm.met <- merge(ea.gmm.met, gmm.names[,.(modules,Name)], by="modules",all.x=T)
  # ea.gmm.met[,modules := Name]
  
  # osa.gmm <- unique(ea.gmm.met$modules)
  
  # Import results on the association with SBP/DBP/HbA1c 
  res.mgs.bp <- fread(paste0(results.folder, 'cor_mgs.gmm_bphb.tsv'))
  
  
  # Import enrichment of GMM ####
  # res.gmm = fread(paste0(results.folder,"ea_GMM.tsv"))
  # res.gmm <- res.gmm[exposure=="t90" & q.value<.05,]
  
  # Positive and negative associations metabolites-species associations 
  
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
  #neg.t90 <- temp.data[rho_t90<0 & q.value_t90<.05 & q.value_odi>=.05, MGS]
  #neg.t90.odi <- temp.data[rho_odi<0 & q.value_t90<.05 & q.value_odi<.05, MGS] # zero
  #neg.odi <- temp.data[rho_odi<0 & q.value_t90>=.05 & q.value_odi<.05, MGS]
  
  #neg.mgs <- c(neg.t90,neg.t90.odi,neg.odi)
  
  #mgs.rel <- c(pos.t90,pos.t90.odi,pos.odi,neg.t90,neg.t90.odi,neg.odi)
  
 
  # Enriched pathways/metabolite groups in the POSITIVE correlations ####
  res.table <- ea_gutsy[Direction=="Positive", .(MGS, Metabolite_subclass, Estimate, q_value)]
  
  # append with gmm enrichment results ####
  # ea.gmm.met[,subclass := fix.subclass.name.fun(subclass)]
  # setnames(ea.gmm.met, c("modules", "subclass", "estimate", "q.value"),
    #       c("MGS","Metabolite_subclass", "Estimate", "q_value"))
  
  # res.table <- rbind(res.table, ea.gmm.met[direction=="positive",.(MGS, Metabolite_subclass, Estimate, q_value)])
  
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
 
  res.table <- ea_gutsy[Direction=="Negative", .(MGS, Metabolite_subclass, Estimate, q_value)]
  
  # append with gmm results 
  # res.table <- rbind(res.table, ea.gmm.met[direction=="negative",.(MGS, Metabolite_subclass, Estimate, q_value)])

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

    
    # Association with T90 and ODI 
    
    res.bm[,exposure:= factor(exposure, levels = c("t90","odi","ahi"))]
    
     t90.odi.cor <- res.bm[MGS %in% pos.mgs & exposure %in% c("t90","odi"), 
                           .(MGS,exposure,rho)]
     t90.odi.q <- res.bm[MGS %in% pos.mgs & exposure %in% c("t90","odi"), 
                           .(MGS,exposure,q.value)]
     
     t90.odi.cor[,exposure := factor(toupper(exposure),levels=c("T90","ODI"))]
     t90.odi.q[,exposure := factor(toupper(exposure),levels=c("T90","ODI"))]
     
     t90.odi.cor <- t90.odi.cor %>% spread(key=exposure, value = rho)
     t90.odi.q <- t90.odi.q %>% spread(key=exposure, value = q.value)

     # create empty rows for gmm 
     # empty.rows <- data.table(MGS=osa.gmm, T90 = 0 , ODI = 0 )
     # t90.odi.cor <- rbind(t90.odi.cor,empty.rows)
     
     # empty.rows <- data.table(MGS=osa.gmm, T90 = 1 , ODI = 1 )
     # t90.odi.q <- rbind(t90.odi.q,empty.rows)
     
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
    assert.order <- function(x){
      return(x[match(pos.mgs, rownames(x)),])
    }
    
  hm.matrix_pos1 <- assert.order(hm.matrix_pos1)
  hm.matrix_neg1 <- assert.order(hm.matrix_neg1) 
  
  qvalues_pos <- qvalues_pos[match(pos.mgs,rownames(qvalues_pos)),]
  qvalues_neg <- qvalues_neg[match(pos.mgs,rownames(qvalues_neg)),]
  
  t90.odi.cor1 <- t90.odi.cor[match(pos.mgs, rownames(t90.odi.cor)),]
  t90.odi.q1 <- t90.odi.q[match(pos.mgs, rownames(t90.odi.q)),]


  # Fix species names to better visualization in the heatmap
  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(hm.matrix_pos1)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(hm.matrix_pos1) <- noms

  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(hm.matrix_neg1)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(hm.matrix_neg1) <- noms

  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(t90.odi.cor1)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(t90.odi.cor1) <- noms

  
  # Merge hm.matrix_pos1 and hm.matrix_neg1
   
  merged.matrix1 <- cbind(hm.matrix_pos1,(hm.matrix_neg1*-1))
  merged.matrix.q1 <- cbind(qvalues_pos,qvalues_neg)
  
  # Legend ####
  col_fun = circlize::colorRamp2(c(0,1.6,2.5),c("white","indianred2", muted("red3")))
  col_fun2 = circlize::colorRamp2(c(0,1.6,2.5),c("white","dodgerblue1", muted("blue4")))
  col_fun3 = circlize::colorRamp2(c(-0.1, 0, 0.1), c(muted("darkcyan"),"white","goldenrod2"))
  
  
  lgd = Legend(col_fun = col_fun, title = "NES", at = c(0,1,2,3), direction = "horizontal",
               labels_gp = gpar(fontsize=6))
  lgd2 = Legend(col_fun = col_fun2, at = c(0,1,2,3), direction = "horizontal",
                labels_gp = gpar(fontsize=6))
  lgd3 = Legend(col_fun = col_fun3, at = c(-.1,-0.05,0,0.05,0.1), direction = "horizontal", 
                title = expression(rho), legend_height = unit(20,"mm"), labels_gp = gpar(fontsize=7))
  
  pd = packLegend(lgd, lgd2, direction = "vertical")
  
  # Annotations ####
  
  ha11 =  HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white", col = 0),
                                            labels = c("negative", "positive"),
                                            labels_gp = gpar(col = "black", fontsize = 8)))
  


      # Create heatmap ####
  set.seed(10)
  ht_opt$TITLE_PADDING = unit(c(.5, .5), "mm")
  ht_opt$COLUMN_ANNO_PADDING = unit(.2, "mm")
  ht_opt$DIMNAME_PADDING = unit(.2, "mm")
  h1 =  Heatmap(merged.matrix1, 
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(merged.matrix.q1[i, j] < 0.05) {
                 grid.text('*', x, y)
               } },
             column_names_rot =  40, 
             column_names_side = "bottom",
             cluster_columns = T,
             cluster_rows = T,   
             show_column_dend = F,
             show_row_dend = FALSE,
             
             column_split = factor(c(rep("positive", ncol(hm.matrix_pos1)),
                                    rep("negative",ncol(hm.matrix_neg1))), 
                                   levels = c("positive","negative")),
             column_gap = unit(2, "mm"),
             top_annotation = ha11,

             #row_title = c("Positive T90/ODI-associated species",
                         #  "Negative T90/ODI-associated species"),
             row_title_gp = gpar(fontsize=10, fill="white",border="black"),
             row_gap = unit(3, "mm"),
             
             show_row_names = FALSE,
             #row_names_side = "left",

             col = circlize::colorRamp2(c(-2.5,-1.6,0,1.6,2.5),c(muted("blue4"), "dodgerblue1", "white",
                                                             "indianred2",muted("red3"))),
             # name="NES (pos)",
             column_labels = colnames(merged.matrix1),
             column_names_gp = gpar(fontsize = 8),
             row_names_gp = gpar(fontsize = 8),
             column_title = "Metabolites associations",
             column_title_gp = gpar(fontsize = 10,fontface='bold'),
             #heatmap_legend_param = list( title = "NES",
                                        # labels_gp = gpar(fontsize=7),
                                        # title_gp = gpar(fontsize=9, fontface="bold"),
                                        # grid_width = unit(2, "mm"), 
                                        # direction = "horizontal")  ,
             width = unit(5*ncol(merged.matrix1),"mm"),
             height = unit(5*nrow(merged.matrix1),"mm"))  

  
 
 # message("Saving plots")
#  dev.off()
  
  # Draw ####
  pdf(file = paste0(results.folder, "mgs_subpathway_ea_heatmap_gutsy_pos_simp.pdf"), 
      width = 11, height = 9)
  set.seed(7)
  draw(h1, heatmap_legend_side = "right", show_heatmap_legend=F, ht_gap = unit(c(3, 6), "mm"))
  draw(pd, x = unit(.96, "npc"), y = unit(.5, "npc"), just = c("right"))
  #draw(bar, x= unit(.371, "npc"), y= unit(.826, "npc"), width = unit(.163, "npc"), height = unit(.18, "npc"),
  #     angle=270, gp = gpar(lwd=.01, col="black"))
  dev.off()
  

  message("Plots saved")
  
  message("Save row order")
  
#  set.seed(10)
#  ht_list <- draw(h1)
#  row.order <- row_order(ht_list)
  
#  row.names.order <- list(
#    positive = row.names(hm.matrix_pos1)[row.order[[1]]], 
#    negative = row.names(hm.matrix_neg1)[row.order[[2]]]
#  )
  
#  saveRDS(row.names.order, file=paste0(results.folder,"row_names_order_hm.rds"))


  
  
  p1 = grid.grabExpr(draw(h1, show_heatmap_legend=F))
plot_grid(p1)  

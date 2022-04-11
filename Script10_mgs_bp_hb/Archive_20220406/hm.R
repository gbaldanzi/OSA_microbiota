# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# This script wil create the Heatmap showing the relation of the T90/ODI-associated
# species and T90-enriched GMMs  with SBP/DBP/Hb1Ac

# Last update: 2022-03-16

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(tidyverse)

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


  # Import result from the association with SBP/DBP/HbA1c

  res.no.bmi <-  fread(paste0(results.folder, "cor.sig.mgs.gmm_bphb.R"))
  res.bmi <- fread(paste0(results.folder, "cor.bmi.sig.mgs.gmm_bphb.R"))
  
  # Separate results from association with species from associations with GMMs
  
  res.mgs.no.bmi <- res.no.bmi[!grep("MF0", MGS_features)]
  res.gmm.no.bmi <- res.no.bmi[grep("MF0", MGS_features)]
  
  res.mgs.bmi <- res.bmi[!grep("MF0", MGS_features)]
  res.gmm.bmi <- res.bmi[grep("MF0", MGS_features)]
  
  # Create matrixes ####
  
  hm.matrix.mgs.no.bmi <- res.mgs.no.bmi[,.(MGS_features, outcomes, rho)] %>% 
    spread(key=outcomes, value = rho)
  setnames(hm.matrix.mgs.no.bmi, c("MGS_features", "DBP_Mean", "Hba1cFormattedResult", "SBP_Mean"), 
           c("MGS", "DBP", "HbA1c", "SBP"))
  setcolorder(hm.matrix.mgs.no.bmi, neworder = c("MGS", "SBP", "DBP", "HbA1c"))
  
  hm.matrix.mgs.bmi <- res.mgs.bmi[,.(MGS_features, outcomes, rho)] %>% 
    spread(key=outcomes, value = rho)
  setnames(hm.matrix.mgs.bmi, c("MGS_features", "DBP_Mean", "Hba1cFormattedResult", "SBP_Mean"), 
           c("MGS", "DBP_BMI", "HbA1c_BMI", "SBP_BMI"))
  setcolorder(hm.matrix.mgs.bmi, neworder = c("MGS", "SBP_BMI", "DBP_BMI", "HbA1c_BMI"))
  
  hm.matrix.mgs <- merge(hm.matrix.mgs.no.bmi,hm.matrix.mgs.bmi, by="MGS")
  
  setDF(hm.matrix.mgs)
  rownames(hm.matrix.mgs) <- hm.matrix.mgs$MGS
  hm.matrix.mgs$MGS <- NULL
  hm.matrix.mgs <- as.matrix(hm.matrix.mgs)
  
  # q-values ####
  q.value.mgs.no.bmi <- res.mgs.no.bmi[,.(MGS_features, outcomes, q.value)] %>% 
    spread(key=outcomes, value = q.value)
  setnames(q.value.mgs.no.bmi, c("MGS_features","DBP_Mean", "Hba1cFormattedResult", "SBP_Mean"), 
           c("MGS", "DBP", "HbA1c", "SBP"))
  setcolorder(q.value.mgs.no.bmi, neworder = c("MGS", "SBP", "DBP", "HbA1c"))
  
  q.value.mgs.bmi <- res.mgs.bmi[,.(MGS_features, outcomes, q.value)] %>% 
    spread(key=outcomes, value = q.value)
  setnames(q.value.mgs.bmi, c("MGS_features", "DBP_Mean", "Hba1cFormattedResult", "SBP_Mean"), 
           c("MGS", "DBP_BMI", "HbA1c_BMI", "SBP_BMI"))
  setcolorder(q.value.mgs.bmi, neworder = c("MGS", "SBP_BMI", "DBP_BMI", "HbA1c_BMI"))
  
  q.value.mgs <- merge(q.value.mgs.no.bmi,q.value.mgs.bmi, by="MGS")
  
  setDF(q.value.mgs)
  rownames(q.value.mgs) <- q.value.mgs$MGS
  q.value.mgs$MGS <- NULL
  q.value.mgs <- as.matrix(q.value.mgs)
  
  # Fix row names 
  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(hm.matrix.mgs)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(hm.matrix.mgs) <- noms
  
  noms <- as.data.table(do.call(rbind, strsplit( split = "____" , rownames(q.value.mgs)) ) )
  noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
  noms <- fix.species.name.fun(noms$V1) 
  rownames(q.value.mgs) <- noms
  
  
  # Assert that rownames are in the proper order (same as previous Heatmap) ####
  
  # Import roworder from previous heatmap
  row.order.hm.list <- readRDS(paste0(results.folder,"row_names_order_hm.rds"))
  row.order.hm <- do.call(c, row.order.hm.list)
  
  hm.matrix.mgs <- hm.matrix.mgs[match(row.order.hm, rownames(hm.matrix.mgs)),]
  q.value.mgs <- q.value.mgs[match(row.order.hm, rownames(q.value.mgs)),]
  
  
  
  set.seed(10)
  ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
  h1 = Heatmap(hm.matrix.mgs, 
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(q.value.mgs[i, j] < 0.05) {
                   grid.text('*', x, y)
                 } },
               column_names_rot=0, 
               column_names_side = "top",
               cluster_columns = FALSE, 
               show_column_dend = F, 
               column_names_gp = gpar(fontsize = 10),
               
               column_split = factor(c(rep("No BMI adjustment",3),
                                    rep("BMI adjusted",3), 
                                  levels = c("No BMI adjustment", "BMI adjusted"))),
               column_gap = unit(2, "mm"),
               
               row_split = factor(c(rep("Positive T90/ODI-associated species",length(row.order.hm.list$positive)),
                                    rep("Negative T90/ODI-associated species",length(row.order.hm.list$negative))), 
                                  levels = c("Positive T90/ODI-associated species", 
                                             "Negative T90/ODI-associated species")),
               cluster_row_slices = F,
               cluster_rows = F,   
               
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
               width = unit(5*ncol(hm.matrix.mgs),"mm")) 
  
  pdf(file = paste0(wrf, "hm.test.pdf"), 
      width = 8, height = 9)
  set.seed(10)
  draw(h1, heatmap_legend_side = "right", merge_legend=F, column_title = "Species-metabolites associations:")
  dev.off()
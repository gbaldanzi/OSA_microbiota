# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Last Update: - 2022-02-25

# Co-abundance network based on the Spearman-correlation between 
# species associated with AHI, T90 or ODI 

# Load packages 
  library(data.table)
  library(corrr)
  library(tidyverse)
  library(ggraph)
  library(igraph)
  
  fix.species.name.fun <- function(char){
    noms <- as.data.table(do.call(rbind, strsplit( split = "____" , char) ) )
    noms[grep("_sp.",V1),V1:=paste0(V1, " (", V2, ")" )]
    noms[grep("_obeum",V1),V1:=paste0(V1, " (", V2, ")" )]
    char2 <- noms$V1
    char2 <- gsub("AF22_8LB","AF22-8LB",char2)
    char2 <- gsub("AM42_11","AM42-11",char2)
    char2 <- gsub("TF06_15AC","TF06-15AC",char2)
    char2 <- gsub("AF46_10NS","AF46-10NS",char2)
    char2 <- gsub("AF36_15AT","AF36-15AT",char2)
    char2 <- gsub("e_P4005","e-P4005", char2)
    char2 <- gsub('[(]HG3A.0168[)]',"",char2)
    char2 <- gsub("_"," ",char2)
    return(char2)
  }
  

# Folders
  results.folder = "/home/baldanzi/Sleep_apnea/Results/"
  input = "/home/baldanzi/Datasets/sleep_SCAPIS/"

# Import data
  pheno <- readRDS(paste0(input,"pheno.MGS.Upp.rds"))
  
  # Import results 
  res.nobmi <- fread(paste0(results.folder, "cor_all.var_mgs.tsv")) 
  res.BMI <- fread(paste0(results.folder, "cor.bmi_all.var_mgs.tsv")) 
  
  mgs.nobmi <- unique(res.nobmi[q.value<.05,MGS])
  mgs.BMI <- unique(res.BMI[q.value<.05,MGS])

  mgs.BMI.pos <- unique(res.BMI[q.value<.05 & rho>0,MGS])
  mgs.BMI.neg <- unique(res.BMI[q.value<.05 & rho<0,MGS])
  
  
  # Correlation between species 
  cor.matrix <- pheno %>% select(all_of(mgs.nobmi)) %>% 
    correlate() %>% stretch()
  
  graph_cors <- cor.matrix %>%
    filter(abs(r) > .1) %>% 
    filter(x %in% mgs.BMI) %>% 
    filter(y %in% mgs.BMI)
  
  col.vector = NULL
  for(i in 1:length(unique(graph_cors$x))){
    col.vector[i] = if(unique(graph_cors$x)[i] %in% mgs.BMI.pos){
      "red"} else if(unique(graph_cors$x)[i] %in% mgs.BMI.neg){
        "blue"}else{"black"}
  }
  
  graph_cors$x <- fix.species.name.fun(graph_cors$x)
  graph_cors$y <- fix.species.name.fun(graph_cors$y)
  
  
  igraph <- graph_cors %>%
    graph_from_data_frame(directed = FALSE)
  

  # Plot
  p = ggraph(igraph) +
    geom_edge_link(aes(edge_alpha = abs(r), 
                       color = r)  ) +
    geom_node_point(color = col.vector) +
    geom_node_text(aes(label = name), 
                   repel = T, size=2) +
    theme_graph() + 
    guides(edge_alpha = "none", edge_width = "none") +
    scale_edge_colour_gradientn(colours = c("dodgerblue2","gray90","firebrick2"),
                                na.value="white", 
                                limits = c(-.2,.4),
                                breaks = c(-.2,0,.2,.4), 
                                values = c(0,0.35,1))
                        

  ggsave(plot = p , filename = "test.png")
# Script 9 - Enrichment analysis plot

# Gabriel Baldanzi 

# Last update: 2022-02-24

  # Packages 
  library(data.table)
  library(tidyverse)
  library(stringr)

  # Folders 
  results.folder = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = paste0(results.folder,"Plots/")
  
  # Import results 
  res.pos = fread(paste0(results.folder,"ea_GMM_pos_new.tsv"))
  
  # Restrict to T90 associations 
  res.pos <- res.pos[exposure=="t90" & q.value<.05,]
  
  # Import pathways/modules names 
  gmm.names <- fread('/home/baldanzi/Datasets/MGS/GMM_reference.csv')
  
  gmm.names[,Name:=str_to_title(Name)]
  gmm.names[,Name:=gsub("Ii","II",Name)]
  
  res.pos <- merge(res.pos,gmm.names,by.x="pathway",by.y="Module",all.x=T,all.y=F)
  
  # Bar plot horizontal 
  res.pos[, Name := factor(Name, levels = Name[order(NES)])]
  
  p <- ggplot(res.pos, aes(y=Name, x=NES)) + 
    geom_bar(aes(fill=q.value), stat="identity", width = .6) + 
    xlab("Normalized enrichment score") + ylab("") +
    scale_x_continuous(expand=c(0,0)) +
    #facet_grid(rows = vars(HL1), scales = "free", switch="y",
     #          space = "free_y") +
    #coord_flip() + 
    theme_classic() + 
    theme(strip.placement = "outside") + 
    scale_fill_continuous(name = expression("p-value"["adj"]))
   # guides(fill = guide_legend())
    
  
  ggsave(plot=p, filename = "GMM_barplot.png", path=output.plot)
  
  ggsave(plot = p, filename = "GMM_barplot.pdf")
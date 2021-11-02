# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 
# Script created in 2021-10-21

# Inferential Statistics 

# This code will produce a Venn Diagram from the correlations of AHI, BMI, 
# and T90 with MGSs in the sensitivity analysis excluding antibiotic users 
# in the last 3 months 

  # Loading packages 
  pacman::p_load(ggplot2,ggvenn, tidyr)
  
  # input and output folders 
  input = "/home/baldanzi/Sleep_apnea/Results/"
  output.plot = "/home/baldanzi/Sleep_apnea/Results/Plots/"

  # Importing results 
  res <- fread(paste0(input,"cor_sa_atb3m_all.var_mgs.tsv"))
  
  res.list <- list(AHI = res[exposure=="ahi",],
                   T90 = res[exposure=="t90",],
                   BMI = res[exposure=="BMI",])

  
  # filter MGS significant at the FDR p-value<0.05
  res.list <- lapply(res.list, function(x){x[q.value>=0.001, q.value:=round(q.value, digits = 3)]})
  
  mgs.fdr = lapply(res.list,function(x){x[x$q.value<0.05,MGS]})
  
  # Create VennDiagram    
  venn <- ggvenn(mgs.fdr,
                    fill_color = c("orange","cornflowerblue","green4"),
                    stroke_color = "white",
                    stroke_size = .2, 
                    show_percentage = F, 
                    fill_alpha = .6) +
    ggtitle("MGS correlated with AHI, T90%, and BMI\nexcluded antibiotic use last 3m") +
    theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5))
  
  ggsave("Venn_bmiahi_sa_atb.png",plot = venn,device = "png", path=output.plot)


  ggsave("Venn_bmiahi_sa_atb.png",plot = venn,device = "png", 
         path="/castor/project/proj_nobackup/wharf/baldanzi/baldanzi-sens2019512")

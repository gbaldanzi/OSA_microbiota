# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Inferential Statistics 

# Version 1: 2021-11-01
# Last Update: - 2021-11-01

# This code will investigate the beta-diversity (Bray Curtis Dissimilarity) 
# in relation to OSA severity groups in 3 different models. 
# After running the PERMANOVA test, this code will now examine the homogeneity of distance
#variance 


# Loading packages 
pacman::p_load(data.table, vegan)


output = "/home/baldanzi/Sleep_apnea/Results/"


# Importing data
pheno <- readRDS("/home/baldanzi/Datasets/sleep_SCAPIS/pheno.MGS.Upp.rds")

pheno <- pheno[valid.ahi=="yes",]

# Transforming two-level factor variables into numeric variables 
a= c("Sex","ppi","metformin","hypermed","dyslipmed")
pheno[,(a):=as.data.frame(data.matrix(data.frame(unclass(pheno[,a, with=F]))))]

# Create dataset containing exclusively the MGS as columns/variables 
mgs.names=grep("____",names(pheno),value=T)
MGSmatrix=pheno[valid.ahi=='yes',mgs.names, with=F]
rownames(MGSmatrix) <- pheno[valid.ahi=="yes", SCAPISid]

# Estimating the BC index - creates a matrix with the BC index between individual samples 
# MGS as relative abundances
BCdist=vegdist(MGSmatrix,method="bray")


#Betadisper
multivariate.dispersions <- betadisper(BCdist,group=pheno[valid.ahi=='yes',OSAcat]) 

saveRDS(multivariate.dispersions,file=paste0(output,"multivariate.dispersions.res.rds"))
# multivariate.dispersions <- readRDS(paste0(output,"multivariate.dispersions.res.rds"))

# Calculate difference in dispersion 
  res2 <- permutest(multivariate.dispersions, pairwise=F, permutations=9999)

  res <- anova(multivariate.dispersions)

  dispersion.df <- data.frame(distance = multivariate.dispersions$distance, 
                              group = multivariate.dispersions$group)
  
  list.res = vector(mode = "list",length = 0)
  
  g1 = "no OSA" ; g2 = c("Mild","Moderate","Severe")
  
  for(i in 1:3){
    message(paste("no_OSA",g2[i],sep = "_"))
    list.res[[paste("no_OSA",g2[i],sep = "_")]] <-  aov(distance ~ group, 
                                                  data = dispersion.df[dispersion.df$group %in% c(g1,g2[i]),])
  }
  
  g1 = "Mild" ; g2 = c("Moderate","Severe")
  
  for(i in 1:2){
    message(paste(g1,g2[i],sep = "_"))
    list.res[[paste(g1,g2[i],sep = "_")]] <-  aov(distance ~ group, 
                                                  data = dispersion.df[dispersion.df$group %in% c(g1,g2[i]),])
  }
  
  g1 = "Moderate" ; g2 = c("Severe")
  
  for(i in 1:1){
    message(paste(g1,g2[i],sep = "_"))
    list.res[[paste(g1,g2[i],sep = "_")]] <-  aov(distance ~ group, 
                                                  data = dispersion.df[dispersion.df$group %in% c(g1,g2[i]),])
  }
  

    clean.res <- function(x){
      res <- lapply(x,function(xx){
      
      p.value <- anova(xx)[[1, "Pr(>F)"]]
      F.value <- anova(xx)[[1, "F value"]]
      df <- data.frame(Comparison = vector("character", length=1),
                       F.value = F.value, p.value = p.value)
      return(df)
      
          })
      res <- do.call(rbind,res)
      res$Comparison <- names(x)
      res$q.value <- p.adjust(res$p.value, method="BH")
      return(res)
    }
 
    
    table.res <-  clean.res(list.res)
    
    fwrite(table.res, paste0(output, "pairwise.permadisp_osa_bc.tsv"))
    




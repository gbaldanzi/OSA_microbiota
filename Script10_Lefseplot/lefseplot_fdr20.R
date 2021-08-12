# # Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Plot results 

# This script will create a Lefse plot based on MGS correlated in model 2 to AHI , T90 or BMI 

  rm(list=ls())
  
  output.plot = "/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/"

  library(data.table)

  #Taxonomy 
  taxonomy = fread("/home/baldanzi/Datasets/MGS/taxonomy", data.table=F)
  
  # Uploading MGS data and restricting to Uppsala participants 
  data.MGS = fread("/home/baldanzi/Datasets/MGS/clean/MGS_relative_abundance_4839_upp_4980_malmo.tsv",data.table=F)
  data.MGS = data.MGS[data.MGS$SCAPISid %in% grep("5-",data.MGS$SCAPISid, value = T ),] 
  
  noms=data.MGS$sample.id
  data.MGS=data.MGS[,grep("____HG3A",names(data.MGS),value=T)]
  rownames(data.MGS)=noms
  
  data.MGS=data.frame(t(data.MGS),stringsAsFactors = F)
  data.MGS$maintax_mgs=rownames(data.MGS)
  
  data=merge(taxonomy,data.MGS,by="maintax_mgs")
  
  # Results from the Spearman correlations 
  
  input <- '/home/baldanzi/Sleep_apnea/Results/'
  
  res.ahi=fread(paste0(input,"cor2_ahi_mgs.tsv"))
  res.t90=fread(paste0(input,"cor2_t90_mgs.tsv"))
  res.bmi=fread(paste0(input,"cor2_BMI_mgs.tsv"))
  
  res = rbind(res.ahi,res.t90,res.bmi)
  
  sig.mgs = unique(res[p.value<.05,][["MGS"]])

  # Restricting annotation to relevant MGS 
  data = data[data$maintax_mgs %in% sig.mgs,]

  data2=data
  
  data$ID=paste("k__",data$superkingdom,"|p__",data$phylum,"|c__",data$class,"|o__",data$order,
                "|f__",data$family,"|g__",data$genus,"|s__",data$maintax_mgs,sep="")
  
  data=data[,c(ncol(data),1:(ncol(data)-1))]
  data=data[,c(1,14:ncol(data))]
  
  data2[data2$family%in%"Oscillospiraceae","family"]="Ruminococcaceae"
  

  mgs.collapse<-function(dades,domain){
    
    variable=names(table(dades[,domain]))
    dades.x=lapply(variable,function (x)   colSums(dades[dades[,domain]%in%x,grep("gut_",names(dades),value=T)]))
    dades.x=data.frame(t(do.call(rbind,dades.x)),stringsAsFactors = F)
    names(dades.x)=variable
    # names(dades.x)=tolower(variable)
    return(dades.x)
    
  }
  
  sup=mgs.collapse(data2,"superkingdom")
  phyl=mgs.collapse(data2,"phylum")
  clas=mgs.collapse(data2,"class")
  ord=mgs.collapse(data2,"order")
  fam=mgs.collapse(data2,"family")
  gen=mgs.collapse(data2,"genus")
  sp=mgs.collapse(data2,"species")
  
  
  sup=sup[names(data)[2:ncol(data)],]
  sup=as.data.frame(t(sup))
  #names(sup)=names(data)[2:ncol(data)]
  phyl=phyl[names(data)[2:ncol(data)],]
  clas=clas[names(data)[2:ncol(data)],]
  ord=ord[names(data)[2:ncol(data)],]
  fam=fam[names(data)[2:ncol(data)],]
  gen=gen[names(data)[2:ncol(data)],]
  sp=sp[names(data)[2:ncol(data)],]
  
  names(sup)=gsub(" ",".",names(sup))
  names(phyl)=gsub(" ",".",names(phyl))
  names(clas)=gsub(" ",".",names(clas))
  names(ord)=gsub(" ",".",names(ord))
  names(fam)=gsub(" ",".",names(fam))
  names(gen)=gsub(" ",".",names(gen))
  names(sp)=gsub(" ",".",names(sp))
  
  
  
  data[nrow(data)+1,]=NA
  data[nrow(data),1]="k__Bacteria"
  data[nrow(data),2:ncol(data)]=sup[1,]
  
  for(i in names(phyl))
  {
    
    n=which(gsub(" ",".",data2$phylum)%in%i)
    for(j in n)
    {
      data[nrow(data)+1,]=NA
      data[nrow(data),1]=paste("k__",gsub(" ",".",data2$sup[j]),"|p__",gsub(" ",".",data2$phyl[j]),sep="")
      data[nrow(data),2:ncol(data)]=phyl[,i]
      
    }
  }
  
  for(i in names(clas))
  {
    n=which(gsub(" ",".",data2$class)%in%i)
    for(j in n)
    {
      data[nrow(data)+1,]=NA
      data[nrow(data),1]=paste("k__",gsub(" ",".",data2$sup[j]),"|p__",gsub(" ",".",data2$phyl[j]),"|c__",gsub(" ",".",data2$class[j]),sep="")
      data[nrow(data),2:ncol(data)]=clas[,i]
      
    }                     
  }
  
  for(i in names(ord))
  {
    n=which(gsub(" ",".",data2$order)%in%i)
    for(j in n)
    {
      data[nrow(data)+1,]=NA
      data[nrow(data),1]=paste("k__",gsub(" ",".",data2$sup[j]),"|p__",gsub(" ",".",data2$phyl[j]),"|c__",gsub(" ",".",data2$class[j]),"|o__",gsub(" ",".",data2$order[j]),sep="")
      data[nrow(data),2:ncol(data)]=ord[,i]
      
    }                     
  }
  
  for(i in names(fam))
  {
    n=which(gsub(" ",".",data2$family)%in%i)
    for(j in n)
    {
      data[nrow(data)+1,]=NA
      data[nrow(data),1]=paste("k__",gsub(" ",".",data2$sup[j]),"|p__",gsub(" ",".",data2$phyl[j]),"|c__",gsub(" ",".",data2$class[j]),"|o__",gsub(" ",".",data2$order[j]),"|f__",gsub(" ",".",data2$family[j]),sep="")
      data[nrow(data),2:ncol(data)]=fam[,i]
      
    }
    
  }
  
  for(i in names(gen))
  {
    n=which(gsub(" ",".",data2$genus)%in%i)
    for(j in n)
    {
      data[nrow(data)+1,]=NA
      data[nrow(data),1]=paste("k__",gsub(" ",".",data2$sup[j]),"|p__",gsub(" ",".",data2$phyl[j]),"|c__",gsub(" ",".",data2$class[j]),"|o__",gsub(" ",".",data2$order[j]),"|f__",gsub(" ",".",data2$family[j]),"|g__",gsub(" ",".",data2$genus[j]),sep="")
      data[nrow(data),2:ncol(data)]=gen[,i]
    }
  }
  
  for(i in names(sp))
  {
    n=which(gsub(" ",".",data2$species)%in%i)
    for(j in n)
    {
      data[nrow(data)+1,]=NA
      data[nrow(data),1]=paste("k__",gsub(" ",".",data2$sup[j]),"|p__",gsub(" ",".",data2$phyl[j]),"|c__",gsub(" ",".",data2$class[j]),"|o__",gsub(" ",".",data2$order[j]),"|f__",gsub(" ",".",data2$family[j]),"|g__",gsub(" ",".",data2$genus[j]),"|s__",gsub(" ",".",data2$species[j]),sep="")
      data[nrow(data),2:ncol(data)]=sp[,i]
      
    }
  }
  
  dat2=unique(data)
  dat=dat2
  
 #  install.packages("/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/microbiomeViz-master",repos = NULL, type = "source")
  library(microbiomeViz)
  
  dat <- data.frame(V1=dat[,1], V2=rowMeans(dat[,-1]), stringsAsFactors = FALSE)
  
  dat=unique(dat)
  dat=dat[order(dat$V1),]
  
  # Create key: MGS and Taxonomic annotation 
  
  splitted <- strsplit(dat2[,1],"[|]")
  
  lastname <- lapply(splitted, function(x){
    a <- x[[length(x)]]
  })
  
  HG3A <- lapply(lastname,function(x){
    substr(x,4,nchar(x))
  })
  
  alltaxa <- data.table(fullname = dat2[,1], 
                        lastname = unlist(lastname),
                        HG3A = unlist(HG3A))
  
  # Reducing the Lefse plot only to fdr20 significant MGS
  fdr.ahi <- res.ahi[q.value<0.20,][["MGS"]]
  fdr.t90 <- res.t90[q.value<0.20,][["MGS"]]
  fdr.bmi <- res.bmi[q.value<0.20,][["MGS"]]
  
  mgs.fdr20 <-  c(fdr.ahi, fdr.t90,fdr.bmi)
  mgs.fdr20 <- unique(mgs.fdr20)
  
  fullname_fdr20 = alltaxa[HG3A %in% mgs.fdr20,][["fullname"]]
  
  dat <- dat[dat$V1 %in% fullname_fdr20,]
  
  # Create Lefse plot backbone 
  
  tr <- parseMetaphlanTSV(dat, node.size.offset=2, node.size.scale=0.1)
  
  print("Creating the tree.backbone - this may take a while")
  p <- tree.backbone(tr, size=0.5)
  
  pdf("/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/tree_fdr20.pdf")
  p
  dev.off()
  
  saveRDS(p, file = paste0(output.plot,"treebackbone_fdr20.rds"))
          
          
  # Create annotation for Lefse plot 
  
  fdr.ahi <- res.ahi[q.value<0.05,][["MGS"]]
  fdr.t90 <- res.t90[q.value<0.05,][["MGS"]]
  fdr.bmi <- res.bmi[q.value<0.05,][["MGS"]]
  
  lefse_lists = data.frame(node=c(alltaxa[HG3A %in% fdr.ahi,lastname],
                                  alltaxa[HG3A %in% fdr.t90,lastname],
                                  alltaxa[HG3A %in% fdr.bmi,lastname]),
                           color=c(rep('indianred3',length(fdr.ahi)),
                                   rep('blue',length(fdr.t90)),
                                   rep('violet',length(fdr.bmi))),
                           stringAsFactor=F)
  
  library(dplyr)
  library(ggtree)
  library(ggplot2)
  
  p2 <- clade.anno2(p, lefse_lists, alpha=0.3)
  
  
  pdf("/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/lefseplot_fdr20_step2.pdf")
  p2
  dev.off()
  
  
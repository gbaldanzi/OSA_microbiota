# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 

# Descriptive Statistics - Taxonomy

# Last update: 2021-10-11

# We found some MGS that were annotated as unclassified all the way to the kingdom level. After 
# discussion with CM, we found that they belonged to the order Clostridiales. Such order does not 
# exist anymore. Therefore, their annotation will be updated to Order: unclassified, Class: Clostridia.
# Phylum: Firmicutes. 

  #### Taxonomy ####
  message("")
  message("Taxonomy")
  message("")
  
# Describe the taxonomic composition by different OSA severity groups. 

# Importing data with taxonomic information for every MGS 
taxonomy = fread("/home/baldanzi/Datasets/MGS/clean/taxonomy",data.table = F)
valid.ahi <- pheno[valid.ahi=='yes',]

# Transforming taxonomic levels into factor variables
a = c("species", "genus", "family", 'order', 'class', 'phylum')
for(i in a){
  taxonomy[[i]]=as.factor(taxonomy[[i]])
}

# Subsetting the data to create a data.frame containing only SCAPISid, ahi, OSAcat
#(OSA severity category), and the relative abundance by phylum for each participants 
a= c("SCAPISid","ahi","OSAcat",levels(taxonomy$phylum))
phylum_abundance=matrix(nrow = nrow(valid.ahi), ncol = length(a))
colnames(phylum_abundance)=a
phylum_abundance=as.data.frame(phylum_abundance)

a=c("SCAPISid","ahi","OSAcat")
phylum_abundance[,a]=as.data.frame(valid.ahi[,a,with=F])

for(i in levels(taxonomy$phylum)){ 
  a = taxonomy$maintax_mgs[taxonomy$phylum==i]
  phylum_abundance[[i]]=rowSums(valid.ahi[,a,with=F])
}

# Determining the mean relative abundance of phylum by OSAcat
phylum_mean= phylum_abundance %>% gather(phyla_name,abundance, levels(taxonomy$phylum)) 
phylum_mean= phylum_mean %>% group_by(OSAcat,phyla_name) %>% summarise(average=mean(abundance))
phylum_mean= as.data.frame(phylum_mean)

#Transforming the OSAcat variable into a factor variable (important for the plot)
phylum_mean$OSAcat = factor(phylum_mean$OSAcat, levels=c("no OSA", "Mild", "Moderate", "Severe"))

# There are 14 phyla. To many for a stacked bar plot. 
# Phyla that have an average below the p25 for all categories will be transformed to "others"
a = quantile(phylum_mean$average,c(.25))
for(i in unique(phylum_mean$phyla_name)){
  if(all(phylum_mean$average[phylum_mean$phyla_name==i]<a)){
    phylum_mean$phyla_name[phylum_mean$phyla_name==i]="other"
  }
}
# Creating the bar plot stacked for phyla relative abundance by OSA cat 
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(14)

p1 = phylum_mean %>% ggplot(aes(x=OSAcat, y=average, fill=phyla_name)) +
  geom_bar(position = "stack", stat = "identity", color="black",lwd=.2) +
  ggtitle("Phyla abundance by sleep apnea severity") +
  ylab("Relative abudance") + xlab("Sleep apnea severity") + 
  theme(title=element_text(hjust=.5, size=12),
        axis.title = element_text(size=12)) +
  scale_fill_manual(name= "phyla", values = mycolors)

#Saving the plot
ggsave("phyla_abundance_OSAcat.png", plot = p1, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")

# Subsetting the data to create a data.frame containing only SCAPISid, ahi, OSAcat
#(OSA severity category), and the relative abundance by family for each participants 
a= c("SCAPISid","ahi","OSAcat",levels(taxonomy$family))
family_abundance=matrix(nrow = nrow(valid.ahi), ncol = length(a))
colnames(family_abundance)=a
family_abundance=as.data.frame(family_abundance)

a=c("SCAPISid","ahi","OSAcat")
family_abundance[,a]=as.data.frame(valid.ahi[,a,with=F])

for(i in levels(taxonomy$family)){ 
  a = taxonomy$maintax_mgs[taxonomy$family==i]
  family_abundance[[i]]=rowSums(valid.ahi[,a,with=F])
}

# Determining the mean relative abundance of family by OSAcat
family_mean = family_abundance %>% gather(family_name,abundance, levels(taxonomy$family)) %>%
  group_by(OSAcat,family_name) %>% summarise(average=mean(abundance))
family_mean = as.data.frame(family_mean)


#Transforming the OSAcat variable into a factor variable (important for the plot)
family_mean$OSAcat = factor(family_mean$OSAcat, levels=c("no OSA", "Mild", "Moderate", "Severe"))

# There are 58 phyla. To many for a stacked bar plot. 
# Famlies that have an average below the p60 for all categories will be transformed to "others"
a = quantile(family_mean$average,c(.70))
for(i in unique(family_mean$family_name)){
  if(all(family_mean$average[family_mean$family_name==i]<a)){
    family_mean$family_name[family_mean$family_name==i]="other"
  }
}

# Creating the bar plot stacked for phyla relative abundance by OSA cat 
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(22)

p2 = family_mean %>% ggplot(aes(x=OSAcat, y=average, fill=family_name)) +
  geom_bar(position = "stack", stat = "identity", color="black",lwd=.2) +
  ggtitle("Family abundance by sleep apnea severity") +
  ylab("Relative abudance") + xlab("Sleep apnea severity") + 
  scale_fill_discrete(name = "families") +
  theme(title=element_text(hjust=.5, size=12),
        axis.title = element_text(size=12)) +
  scale_fill_manual(name= "families", values = mycolors)

#Saving the plot
ggsave("family_abundance_OSAcat.png", plot = p2, device = "png", 
       path="/home/baldanzi/Sleep_apnea/Descriptive/")


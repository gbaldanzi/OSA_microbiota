# Project: Sleep apnea and gut microbiota
# Gabriel Baldanzi 


# This script will prepare the data to the Spearman's correlation between 
# GMM abundance and plasma metabolites. 

# This analysis will include participants from both Uppsala and Malmo 

rm(list = ls())

  # load packages 
  library(data.table)
  library(Hmisc)
  library(sjmisc)
  library(dplyr)
  library(tidyr)
  library(vegan)

  # Folders 
  input <- "/home/baldanzi/Datasets/"
  output <- '/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/data_processed/'

  # Uploading phenotype data 
  pheno=fread(paste0(input,"sleep_SCAPIS/SCAPIS-DATA-PETITION-170-20210315.csv"), header = T, na.strings=c("", "NA"))
  
  # Uploading MGS/gut microbiota data 
    data.MGS = fread(paste0(input,'MGS/clean/MGS_relative_abundance_4839_upp_4980_malmo.tsv'))
  
  
    # Fixing MGS names to latest annotation 
    taxonomy = fread(paste0(input,"MGS/taxonomy"))
    mgs.names.index <- grep("____",names(data.MGS))
    names(data.MGS)[mgs.names.index] <- taxonomy$maintax_mgs
   
  # Merge pheno with data.MGS ####
    pheno <- merge(pheno, data.MGS, by="Subject", all=F)
    
  # Clinical microbiomics data ####
    cmvar <- fread("/home/baldanzi/Datasets/MGS/clean/CM_lab_variables_clean.tsv")
    
    a = c("Subject", "plate", "received", "read.pairs.per.sample", "dna.yield")
    pheno = merge(pheno, cmvar[,a,with=F], by="Subject", all.x=T, all.y=F)
    
    # Plate by site 
    pheno[,site_plate := paste(Site,plate,sep="_")]
    
    pheno[,site_plate := as.factor(site_plate)]
    pheno[,received := as.factor(received)]
    pheno <- pheno[,-c("SCAPISid")]
    
    setnames(pheno, "Subject", "SCAPISid")
    
  # Data on coutry/place of birth ####
    pob = fread("/home/baldanzi/Datasets/Placeofbirth.csv", header=T, sep=",")
    pob$placebirth = factor(pob$q005a, levels = c("scandinavia", "europe","asia","other"))
    pheno <- merge(pheno, pob[,.(SCAPISid, placebirth)], by= "SCAPISid", all.x=T, all.y=F)
  
    # Recoding variables ####
    
    # Month of anthropometric collection date - variable to account for seasonality
    pheno[,visit.month:=format(as.POSIXct(pheno$AnthropometryCollectionDate),"%B")]
    pheno[,visit.month:=factor(visit.month,c(month.name,"June.July"))]
    
    # Merging June to July due to the low number of July participants 
    pheno[visit.month %in% c("June","July"),visit.month:="June.July"]
    
    #Age
    pheno[,age:=AgeAtVisitOne]
    
    #Sex
    pheno$Sex = factor(pheno$Sex, levels = c("MALE", "FEMALE"), 
                       labels = c("male", "female"))
    
    #smoking 
    setnames(pheno, 'derived_smoke_status', 'smokestatus')
    pheno$smokestatus = rec(pheno$smokestatus, 
                            rec = "NEVER=never; CURRENT=current; EX_SMOKER=former; NA,UNKNOWN=NA")
    pheno$smokestatus = factor(pheno$smokestatus, levels = c("never", "former", "current"))
    
    # Highest education 
    pheno$educat = rec(pheno$cqed001, rec = "-99=NA; else=copy")
    pheno$educat = factor(pheno$educat, levels = c(0,1,2,3), 
                          labels = c("uncompleted primary or lower secondary", 
                                     "lower secondary education", "upper secondary education", 
                                     "university education"))
    
    # BMI categories 
    pheno$BMIcat = rec(pheno$BMI, rec = "min:24.9=1; 25:29.9=2; 30:max = 3")
    pheno$BMIcat = factor(pheno$BMIcat, levels = c(1,2,3), labels = c("<25","25-30",">=30"))
    
    # Self-reported physical activity 
    pheno$leisurePA = rec(pheno$cqpa012, rec = "4=NA; else=copy")
    pheno$leisurePA = factor(pheno$leisurePA, levels = c(0,1,2,3), 
                             labels = c("mostly sedentary", "moderate activity", 
                                        "regular and moderate activity ", "regular exercise or training"))
    
    # Diabetes variable
    pheno$diabd = rec(pheno$Diabetes, rec = "NORMOGLYCEMIA=0; ELEV_HBA1C,IFG=1; KNOWN_DM,NEW_DM=2")
    pheno$diabd = factor(pheno$diabd, levels = c(0,1,2), 
                         labels = c("normoglycemic", "impaired glucose tolerance", "type 2 diabetes" ))
    
    # Site 
    pheno$Site <-  factor(pheno$Site)
    
    
    # Diet ####
    # Removing over- and under- reported based on the 3SD of the natural log of Energy intake
    
    # remove extremes outliers 
    pheno[Energi_kcal==0, Energi_kcal:=NA]
    pheno[,log.energi:= log(pheno$Energi_kcal)] # Calculating the log energy intake
    
    # Estimate sd and mean of energy intake
    male.sd <-  sd(pheno[Sex=="male",log.energi],na.rm=T)
    male.mean <- mean(pheno[Sex=="male",log.energi],na.rm=T)
    female.sd <- sd(pheno[Sex=="female",log.energi],na.rm=T)
    female.mean <- mean(pheno[Sex=="female",log.energi],na.rm=T)
    
    # Using 3 sd to categorize over- and under-reporters
    pheno[!is.na(Energi_kcal),energi.reporter:="ok"]
    
    ll <- female.mean-(3*female.sd)
    ul <- female.mean+(3*female.sd)
    pheno[Sex=="female" & log.energi<ll, energi.reporter:="under"]
    pheno[Sex=="female" & log.energi>ul, energi.reporter:="over"]
    
    ll <- male.mean-(3*male.sd)
    ul <- male.mean+(3*male.sd)
    pheno[Sex=="male" & log.energi<ll, energi.reporter:="under"]
    pheno[Sex=="male" & log.energi>ul, energi.reporter:="over"]
    
    # Assign NA to diet variables for under- or over- reporter
    pheno[,energi.original:=Energi_kcal]
    pheno[energi.reporter!="ok",Energi_kcal:=NA]
    pheno[energi.reporter!="ok",Fibrer:=NA]
    
    
    # Self-reported hypertension variable
    pheno$hypertension = factor(pheno$cqhe034, levels = c("NO", "YES"), labels = c("no", "yes"))
    
    # Self-reported dyslipidemia 
    pheno$dyslipidemia = factor(pheno$cqhe036, levels = c("NO", "YES"), labels = c("no", "yes"))
    
    # Self-reported diabetes medication 
    pheno$diabmed = factor(pheno$cqhe039, levels = c("NO", "YES"), labels = c("no", "yes"))
    
    # Self-reported hypertension medication 
    pheno$hypermed = factor(pheno$cqhe035,  levels = c("NO", "YES"), labels = c("no", "yes"))
    
    # Self-reported dyslipidemia medication 
    pheno$dyslipmed = factor(pheno$cqhe037,  levels = c("NO", "YES"), labels = c("no", "yes"))
    
    
    # Metabolon data ####
    metabolon=fread("/home/baldanzi/Datasets/Metabolon/clean/scapis_metabolon_batchnorm_scapisid.tsv", header=T)
    a = c("SCAPISid","MET_100002725","MET_100002808","MET_100002405")
    medication = metabolon[,a,with=F]
    medication[,ppi:=ifelse(MET_100002725>min(MET_100002725,na.rm=T) |
                              MET_100002808>min(MET_100002808,na.rm=T),1,0)]
    
    medication[,ppi:=factor(ppi, levels = c(0,1), 
                            labels = c("no", "yes"))]
    
    medication[,metformin:=ifelse(MET_100002405>min(MET_100002405,na.rm=T),1, 0)]
    
    medication[,metformin:=factor(metformin, levels = c(0,1), 
                                  labels = c("no", "yes"))]
    
    a = c("SCAPISid","ppi","metformin")
    pheno <- merge(pheno,  medication[,a,with=F], by="SCAPISid", all.x=T)
    
    
    # GMM results 
    
    #Importing enrichment analysis results 
    results.folder <-  "/proj/nobackup/sens2019512/users/baldanzi/sleepapnea_gut/results/"
    res.ea.gmm <- fread(paste0(results.folder,"ea_GMM.tsv"))
    
    osa.gmm <- unique(res.ea.gmm[q.value<.05,pathway])
    
    # Import GMM abundance ####
    library(rio)
    GMM.uppsala <- import('/proj/sens2019512/SCAPIS_org/SCAPIS/final_release_CMv1/Uppsala/upugut03.GMMComp.percent.xlsx')
    GMM.malmo <- import('/proj/sens2019512/SCAPIS_org/SCAPIS/final_release_CMv1/Malmo/lungut03.GMMComp.percent.xlsx')
    
    GMM.uppsala <- GMM.uppsala[GMM.uppsala$Module %in% osa.gmm, -which(names(GMM.uppsala) %in% "annotations")]
    GMM.malmo <- GMM.malmo[GMM.malmo$Module %in% osa.gmm, -which(names(GMM.malmo) %in% "annotations")]
    
    GMM.table <- merge(GMM.uppsala, GMM.malmo, by=c("Module"))
    names(GMM.table) <- gsub("_b", "", names(GMM.table)) # Making ids compatible with the data.MGS 
    names(GMM.table) <- gsub("b", "", names(GMM.table))
    names(GMM.table) <- gsub("c", "", names(GMM.table))
    
    rownames(GMM.table) <- GMM.table$Module
    GMM.table$Module <- NULL
    
    GMM.table <- t(GMM.table)
    modules.table <- colnames(GMM.table)
    id <- rownames(GMM.table)
    
    GMM.table <- as.data.frame(GMM.table)
    names(GMM.table) <- modules.table
    GMM.table$sample.id <- id 
    
    pheno <- merge(pheno, GMM.table, by="sample.id", all.x=T, all.y=F)
    
    # Import metabolites data 
    input2 <- "/proj/sens2019512/SCAPIS_org/SCAPIS/metabolon_final/clean/"
    metabolites <- fread(paste0(input2,'scapis_merged_data_batchnorm_clean.tsv'))
    annotation <- fread(paste0(input2, 'scapis_merged_annotations_batchnorm_clean.tsv'))
    
    Drug_MET_ID <- annotation[grep("Drug", SUB_PATHWAY) , MET_ID]
    
    Drug_MET_ID <- Drug_MET_ID[-which(Drug_MET_ID=='MET_100020837')]
    
    metabolites[, (Drug_MET_ID) := lapply(.SD, function(x) {
      ifelse(x>min(x,na.rm=T), 1, 0) }), .SDcols = Drug_MET_ID
    ]
    
    pheno[SCAPISid %in% metabolites$SCAPIS_ID, metabolon_data := T]
    
    pheno <- merge(pheno, metabolites, by.x="SCAPISid", by.y="SCAPIS_ID", 
                   all.x=T, all.y=F)
    
    # Calculate shannon index 
    # To keep the Shannon index calculations in the same as performed in the 
    # Gusty Atlas, we will remove the species present in less than 100 individual
    
    species.names <- grep("HG3A",names(pheno), value=T)
    
    non.rare.species <- lapply(species.names, function(x) {
      a <- nrow(pheno[metabolon_data == T & pheno[[x]]!=0,])
      if(a>100){return(x)}
    })
    
    non.rare.species <- do.call(c,non.rare.species)
    
    pheno[, shannon:= diversity(pheno[,non.rare.species,with=F])]
    
    
    # Saving results ####
    saveRDS(pheno, file = paste0(output,"phenotype_upp_malm.rds"))
    
    
      # Prepare subclass list for the enrichment analysis 
    annotation <- fread('/home/baldanzi/Datasets/Metabolon/clean/original/scapis_annotations_clean.tsv')
    annotation[,SUB_PATHWAY:=gsub("-","_",SUB_PATHWAY)]
    annotation[,SUB_PATHWAY:=gsub("/","_",SUB_PATHWAY)]
    annotation[,SUB_PATHWAY:=gsub('(',"_",SUB_PATHWAY,fixed=T)]
    annotation[,SUB_PATHWAY:=gsub(")","_",SUB_PATHWAY)]
    annotation[,SUB_PATHWAY:=gsub(",","",SUB_PATHWAY)]
    annotation[,SUB_PATHWAY:=gsub(";","",SUB_PATHWAY)]
    annotation[,SUB_PATHWAY:=gsub(" ","_",SUB_PATHWAY)]
    lev <- annotation[,SUB_PATHWAY]
    lev <- lev[!is.na(lev)]
    lev <- lev[lev!=""]
    pathways=NULL
    pathways=list()
    for(i in lev)
    {
      eval(parse(text=paste0("pathways$", i , "=annotation[SUB_PATHWAY==i,MET_ID]")))
    }
    
    
    saveRDS(pathways, paste0(output,'subclass.list.rds'))
  
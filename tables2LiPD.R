#load in read in SISAL data and populate LiPD files.
library(here)
library(tidyverse)
library(magrittr)
library(lipdR)
load(here("allTables.Rdata"))

#load in name converter
nc <- read_csv(here("sisal2lipd_names.csv")) %>% 
  filter(!is.na(lipd))

nc.chron <- read_csv(here("sisal2lipd_names_chron.csv")) %>% 
  filter(!is.na(lipd))

nc.chron.lam <- read_csv(here("sisal2lipd_names_chron_lam.csv")) %>% 
  filter(!is.na(lipd))
#assign in data.frames

sites <- all$site#each site
entity <- all$entity #each speleothem/collection?
elr <- all$entity_link_reference #links entities to references
pubs <- all$reference
sample <- all$sample
d18O <- all$d18O
d13C <- all$d13C
dating <- all$dating
datingLamina <- all$dating_lamina
origDates <- all$original_chronology

#loop through sites
for(s in 19:nrow(sites)){
#for(s in 1:10){
  #GEO
  geo <- list()
  geo$siteName <- sites$site_name[s]
  geo$latitude <- sites$latitude[s]
  geo$longitude <- sites$longitude[s]
  geo$elevation <- sites$elevation[s]
  geo$sisalSiteId <- sites$site_id[s]
  geo$geology <- sites$geology[s]
  geo$rockAge <- sites$rock_age[s]
  geo$hasMonitoring <- sites$monitoring[s]
  
  #get entities for site
  this.ent <- dplyr::filter(entity,site_id == geo$sisalSiteId)#connect this site to it's stals
  this.elr <- dplyr::filter(elr,entity_id %in% this.ent$entity_id)#connect these stals to their pubs
  this.pub <- dplyr::filter(pubs,ref_id %in% this.elr$ref_id)#grab the pubs in this site
  
  #count number of unique pubs
  this.elr$thisRefId <- NA
  this.elr$thisRefId[1] <- 1
  if(nrow(this.elr)>1){
    for(i in 2:nrow(this.elr)){
      if(any(this.elr$ref_id[i] %in% this.elr$ref_id[1:(i-1)])){
        this.elr$thisRefId[i] <- this.elr$thisRefId[min(which(this.elr$ref_id[i] == this.elr$ref_id[1:(i-1)]))]
      }else{
        this.elr$thisRefId[i] <- max(this.elr$thisRefId[1:(i-1)])+1
      }
    }
  }
  #PUB
  pub <- vector(mode = "list",length = nrow(this.pub))
  for(p in 1:nrow(this.pub)){
    pub[[p]]$citation <- this.pub$citation[p]
    pub[[p]]$DOI <- this.pub$publication_DOI[p]
    re <- gregexpr(pattern = "[1-2][0-9][0-9][0-9]",pub[[p]]$citation)
    
    pub[[p]]$year <- as.numeric(stringr::str_sub(pub[[p]]$citation,start = re[[1]][1],end = re[[1]][1]+3))
  }
  if(length(pub)>0){#grab the first author for later
    split <- (stringr::str_split(pub[[1]]$citation,","))
    sl <- sapply(split,stringr::str_length)
    firstAuthor <- split[[1]][min(which(sl>=3))]
    pubYear1 <-  pub[[1]]$year
  }
  
  #BASE
  dataSetName <-  stringr::str_remove_all(stringr::str_c(geo$siteName,firstAuthor,pub[[p]]$year,sep = "."),pattern = " ")
  
  #PALEODATA
  pmt <- vector(mode = "list",length = nrow(this.ent))
  for(e in 1:nrow(this.ent)){#create a measurement table for every entry
    pmt[[e]]$tableName <- this.ent$entity_name[e]
    pmt[[e]]$SISALEntityID <- this.ent$entity_id[e]
    pmt[[e]]$hasPublication <- this.elr$thisRefId[e]
    
    #find samples that belong to this entity
    this.samp <- filter(sample,entity_id == this.ent$entity_id[e])
    
    #create a complete data.frame for this entity
    adf <- left_join(this.samp,d18O,by="sample_id") %>% 
      left_join(d13C,by="sample_id") %>% 
      left_join(origDates,by="sample_id") %>% 
      left_join(all$gap,by="sample_id") %>% 
      left_join(all$hiatus,by="sample_id")
    
    
    for(n in 1:nrow(nc)){#create a column for each row
      
      tc <- list()#create an empty list
      #COLUMN META
      tc$variableName <- nc$lipd[n]
      tc$units <- nc$units[n]
      tc$description <- nc$description[n]
      tc$variableType <- nc$variableType[n]
      if(!is.na(nc$proxyObservationType[n])){
        tc$proxyObservationType <- nc$proxyObservationType[n]
      }
      if(!is.na(nc$inferredVariableType[n])){
        tc$inferredVariableType <- nc$inferredVariableType[n]
      }
      
      #add in the data
      tc$values <- as.matrix(adf[nc$sisal[n]])
      
      if(!all(is.na(tc$values))){      #plop into the measuermentTable
        pmt[[e]][[tc$variableName]] <- tc
      }
    }#end column loop
    
  }#end measurementTable loop
  
  #CHRONDATA
  cmt <- vector(mode = "list",length = nrow(this.ent))
  hasChron <- c()
  
  for(e in 1:nrow(this.ent)){#create a measurement table for every entry
    cmt[[e]]$tableName <- this.ent$entity_name[e]
    cmt[[e]]$SISALEntityID <- this.ent$entity_id[e]
    cmt[[e]]$hasPublication <- this.elr$thisRefId[e]
    
    #find dates that belong to this entity
    this.dates <- filter(dating,entity_id == this.ent$entity_id[e])
    
    #find dates that belong to this entity
    this.datingLamina <- filter(datingLamina,entity_id == this.ent$entity_id[e])
    hasChron[e] <- TRUE
    

    if(nrow(this.dates)>0){#tie points
      for(n in 1:nrow(nc.chron)){#create a column for each row
        
        tc <- list()#create an empty list
        #COLUMN META
        tc$variableName <- nc.chron$lipd[n]
        tc$units <- nc.chron$units[n]
        tc$description <- nc.chron$description[n]
        tc$variableType <- nc.chron$variableType[n]
        if(!is.na(nc.chron$proxyObservationType[n])){
          tc$proxyObservationType <- nc.chron$proxyObservationType[n]
        }
        if(!is.na(nc$inferredVariableType[n])){
          tc$inferredVariableType <- nc.chron$inferredVariableType[n]
        }
        
        #add in the data
        col.name <- nc.chron$sisal[n]
        if(grepl("[0-9]",substr(col.name,1,1))){#starts with a number
          #append an X
          col.name <- str_c("X",col.name)
        }
        tc$values <- as.matrix(this.dates[col.name])
        tc$values <- str_remove_all(tc$values,'"')
        tc$values <- str_remove_all(tc$values,"'")
        tc$values <- str_remove_all(tc$values,",")
        
        if(!all(is.na(tc$values))){ #plop into the measuermentTable
          cmt[[e]][[tc$variableName]] <- tc
          tableRows <- length( tc$values)
        }
      }#end column loop
      
      #add an agetype column
      tc <- list()#create an empty list
      #COLUMN META
      tc$variableName <- "ageType"
      tc$units <- "unitless"
      tc$description <- "broad class of age control"
      tc$variableType <- "sampleMetadata"
      
      #add in the data
      tc$values <- matrix("U/Th",nrow =tableRows)
      cmt[[e]]$ageType <- tc
      
    }#end tiepoints
    if(nrow(this.datingLamina)>0){#lamina
      c.lam <- list()
      for(n in 1:nrow(nc.chron.lam)){#create a column for each row
        tc <- list()#create an empty list
        #COLUMN META
        tc$variableName <- nc.chron.lam$lipd[n]
        tc$units <- nc.chron.lam$units[n]
        tc$description <- nc.chron.lam$description[n]
        tc$variableType <- nc.chron.lam$variableType[n]
        if(!is.na(nc.chron.lam$proxyObservationType[n])){
          tc$proxyObservationType <- nc.chron.lam$proxyObservationType[n]
        }
        if(!is.na(nc.chron.lam$inferredVariableType[n])){
          tc$inferredVariableType <- nc.chron.lam$inferredVariableType[n]
        }
        
        #add in the data
        col.name <- nc.chron.lam$sisal[n]
        if(grepl("[0-9]",substr(col.name,1,1))){#starts with a number
          #append an X
          col.name <- str_c("X",col.name)
        }
        tc$values <- as.matrix(this.datingLamina[col.name])
        tc$values <- str_remove_all(tc$values,'"')
        tc$values <- str_remove_all(tc$values,"'")
        tc$values <- str_remove_all(tc$values,",")
        
        if(!all(is.na(tc$values))){ #plop into the measuermentTable
          cmt[[e]][[tc$variableName]] <- tc
          tableRows <- length( tc$values)
        }
      }#end column loop
      
      #add an agetype column
      tc <- list()#create an empty list
      #COLUMN META
      tc$variableName <- "ageType"
      tc$units <- "unitless"
      tc$description <- "broad class of age control"
      tc$variableType <- "sampleMetadata"
      
      #add in the data
      tc$values <- matrix("layerCount",nrow =tableRows)
      c.lam$ageType <- tc
      
    }#end lamina
    
    #combine as needed
    if(nrow(this.dates)>0 & nrow(this.datingLamina)>0){
      tieLength <- length(cmt[[e]]$age$values)
      lamLength <- length(c.lam$age$values)
      tieNames <- c(nc.chron$lipd,"ageType")
      lamNames <- c(nc.chron.lam$lipd,"ageType")
      
      
      for(n in 1:length(tieNames)){#loop through tiepoints and append...
        this.col <- cmt[[e]][[tieNames[n]]]
        if(any(this.col$variableName %in% lamNames)){#append lams
          print(paste("appending",this.col$variableName))
          this.col$values <- c(this.col$values,c.lam[[this.col$variableName]]$values)
        }else{#append NAs
          this.col$values <- c(this.col$values,matrix(NA,nrow=lamLength))
        }
        cmt[[e]][[tieNames[n]]] <- this.col
      }
      for(n in 1:length(lamNames)){#loop through tiepoints and append...
        this.col <- c.lam[[lamNames[n]]]
        if(any(this.col$variableName %in% tieNames)){#append lams
          #do nothing
        }else{#prepend NAs
          this.col$values <- c(matrix(NA,nrow=tieLength),this.col$values)
          cmt[[e]][[lamNames[n]]] <- this.col
        }
      }
    }else if(nrow(this.datingLamina)>0){
      cmt[[e]] <- c.lam
    }
  
    #
    if(nrow(this.dates)==0 & nrow(this.datingLamina)==0){
      hasChron[e] <- FALSE
    }   
  }#end measurementTable loop
  
  
  #build into a lipd file
  L <- list()
  L$dataSetName <- dataSetName
  L$lipdVersion <- 1.3
  L$geo <- geo
  L$pub <- pub
  L$paleoData <- vector(mode = "list",length = 1)
  L$paleoData[[1]]$measurementTable <- pmt
  if(any(hasChron)){
    L$chronData <- vector(mode = "list",length = 1)
    L$chronData[[1]]$measurementTable <- cmt
  }
  
  writeLipd(L,path = here("lipds"),ignore.warnings = TRUE)
  
}

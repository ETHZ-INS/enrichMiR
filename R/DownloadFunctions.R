
# Downlaod Targetscan Position Files
.loadTargetscanPos <- function(species = c("human","mouse","rat")) {
  species <- match.arg(species)
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116, 
                  stop("No matched species"))
  
  # download TargetScan miRNA Positions
  tmp <- tempfile()
  if (species == "human"){
    #Downlaod Targetscan Species specific UTR file
    download.file(
      "http://www.targetscan.org/vert_72/vert_72_data_download/Predicted_Targets_Info.default_predictions.txt.zip", tmp)
    miRPos <- fread(unzip(file.path(tmp)),drop = c("MSA start", "MSA end", "PCT") )
  }else if(any(species %in% c("mouse","rat"))){
    #Downlaod Targetscan Species specific UTR file
    download.file(
      "http://www.targetscan.org/mmu_72/mmu_72_data_download/Conserved_Family_Conserved_Targets_Info.txt.zip", tmp)
    miRPos <- fread(unzip(file.path(tmp)),drop = c("MSA start", "MSA end", "PCT") )
  }
  unlink(tmp)
  
  #filter for the miRNA species
  miRPos <- miRPos[miRPos$`Species ID` == spec,]
  miRPos$`Gene ID` = substr(miRPos$`Gene ID`,1,18)
  miRPos$`Transcript ID` = substr(miRPos$`Transcript ID`,1,18)
  miRPos
}



# Download Targetscan miRNA families
.getMirfamilies <- function(species = c("human","mouse","rat")) {
  species <- match.arg(species)
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116, 
                  stop("No matched species"))
  
  # download TargetScan miRNA families
  tmp <- tempfile()
  download.file("http://www.targetscan.org/vert_72/vert_72_data_download/miR_Family_Info.txt.zip", tmp)
  miRFam <- fread(unzip(file.path(tmp)))
  unlink(tmp)
  
  #filter for the miRNA species
  miRFam <- miRFam[miRFam$`Species ID` == spec,]
  
  
}



# Download Targetscan Conserved Sites
.loadTargetscanSites <- function(species = c("human","mouse","rat")) {
  species <- match.arg(species)
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116, 
                  stop("No matched species"))
  
  # download TargetScan conserved miRNA sites
  tmp <- tempfile()
  if (species == "human"){
    #Downlaod Targetscan Species specific site file
    download.file(
      "http://www.targetscan.org/vert_72/vert_72_data_download/Summary_Counts.default_predictions.txt.zip", tmp)
    ConSites <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
  }else if(any(species %in% c("mouse","rat"))){
    #Downlaod Targetscan Species specific site file
    download.file(
      "http://www.targetscan.org/mmu_72/mmu_72_data_download/Summary_Counts.default_predictions.txt.zip", tmp)
    ConSites <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
  }
  unlink(tmp)
  
  #filter for the miRNA species
  ConSites <- ConSites[ConSites$`Species ID` == spec,]
  ConSites$`Transcript ID` = substr(ConSites$`Transcript ID`,1,15)
  ConSites
}
    
    
    
    

# Download All Targetscan Sites
.loadTargetscanSitesAll <- function(species = c("human","mouse","rat")) {
  species <- match.arg(species)
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116, 
                  stop("No matched species"))
  
  # download TargetScan conserved miRNA sites
  tmp <- tempfile()
  if (species == "human"){
    #Downlaod Targetscan Species specific site file
    download.file(
      "http://www.targetscan.org/vert_72/vert_72_data_download/Summary_Counts.all_predictions.txt.zip", tmp)
    AllSites <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
  }else if(any(species %in% c("mouse","rat"))){
    #Downlaod Targetscan Species specific site file
    download.file(
      "http://www.targetscan.org/mmu_72/mmu_72_data_download/Summary_Counts.all_predictions.txt.zip", tmp)
    AllSites <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
  }
  unlink(tmp)
  
  #filter for the miRNA species
  AllSites <- AllSites[AllSites$`Species ID` == spec,]
  AllSites$`Transcript ID` = substr(AllSites$`Transcript ID`,1,15)
  AllSites
}    
    
    
    
    
    
    
    
    
    
    
    
    
    

#' importFrom Matrix sparseMatrix

.fetch_mirtarbase <- function(species, returnType=c("dataframe","matrix")){
  
  returnType <- match.arg(returnType)
  
  species <- switch(species, human="hsa", mouse="mmu", rat="rno", species)
  
  tmp <- tempfile()
  
  download.file(paste0("http://mirtarbase.cuhk.edu.cn/cache/download/8.0/", 
                       
                       species, "_MTI.xls"), destfile = tmp)
  
  e <- readxl::read_excel(tmp)
  
  e <- e[,c(2,4)]
  
  e <- as.data.frame(e)
  
  e[,1] <- as.factor(e[,1])
  
  e[,2] <- as.factor(e[,2])
  
  if(returnType=="dataframe"){
    
    colnames(e) <- c("set","feature")
    
    return(e)
    
  }
  
  sparseMatrix( as.numeric(e[,2]), as.numeric(e[,1]), 
                
                dimnames=list(levels(e[,2]), levels(e[,1])) )
  
}
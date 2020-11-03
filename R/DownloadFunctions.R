
#' Downlaod Targetscan Position Files
.fetch_TargetscanPos <- function(species = c("human","mouse","rat"),type=c("conserved","all")) {
  species <- match.arg(species)
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116, 
                  stop("No matched species"))
  
  fams <- .fetch_Mirfamilies(species)
  
  
  
  
  
  # # download TargetScan miRNA Positions
  # tmp <- tempfile()
  # if (species == "human"){
  #   #Downlaod Targetscan Species specific UTR file
  #   download.file(
  #     "http://www.targetscan.org/vert_72/vert_72_data_download/Predicted_Targets_Info.default_predictions.txt.zip", tmp)
  #   miRPos <- fread(unzip(file.path(tmp)),drop = c("MSA start", "MSA end", "PCT") )
  # }else if(any(species %in% c("mouse","rat"))){
  #   #Downlaod Targetscan Species specific UTR file
  #   download.file(
  #     "http://www.targetscan.org/mmu_72/mmu_72_data_download/Conserved_Family_Conserved_Targets_Info.txt.zip", tmp)
  #   miRPos <- fread(unzip(file.path(tmp)),drop = c("MSA start", "MSA end", "PCT") )
  # }
  # unlink(tmp)
  # 
  #filter for the miRNA species
  # miRPos <- miRPos[miRPos$`Species ID` == spec,]
  # miRPos$`Gene ID` <- gsub("\\..*","",miRPos$`Gene ID`)
  # miRPos$`Transcript ID` <- gsub("\\..*","",miRPos$`Transcript ID`)
  # miRPos
}



#' Download Targetscan miRNA families
.fetch_Mirfamilies <- function(species = c("human","mouse","rat")) {
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
  miRFam <- read.delim(unzip(tmp), header=TRUE)
  unlink(tmp)
  
  #filter for the miRNA species
  miRFam <- miRFam[miRFam$Species.ID == spec,]
  miRFam
}

.ens2symbol <- function(species=c("human","mouse","rat")){
  species <- match.arg(species)
  symbol <- switch(species,
                   human="hgnc_symbol",
                   mouse="mgi_symbol",
                   rat="rgd_symbol",
                   "uniprot_gn_symbol")
  ensembl <- switch(species,
               human = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl"),
               mouse = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl"),
               rat = useMart("ENSEMBL_MART_ENSEMBL",dataset="rnorvegicus_gene_ensembl"), 
               stop("No matched species"))
  anno <- getBM(attributes=c("ensembl_gene_id",symbol), mart = ensembl)
  anno <- anno[anno[,2]!="",]
  anno <- anno[!duplicated(anno[,1]),]
  e2s <- anno[,2]
  names(e2s) <- anno[,1]
  e2s
}

.getTargetScan <- function(species = c("human","mouse","rat"), type=c("conserved","all"), keepMers=FALSE){
  type <- match.arg(type)
  species <- match.arg(species)
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116, 
                  stop("No matched species"))
  
  fams <- .fetch_Mirfamilies(species)
  
  # download TargetScan conserved miRNA sites
  tmp <- tempfile()
  if(type=="conserved"){
    if (species == "human"){
      #Downlaod Targetscan Species specific site file
      download.file(
        "http://www.targetscan.org/vert_72/vert_72_data_download/Summary_Counts.default_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }else if(any(species %in% c("mouse","rat"))){
      #Downlaod Targetscan Species specific site file
      download.file(
        "http://www.targetscan.org/mmu_72/mmu_72_data_download/Summary_Counts.default_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }
    a$sites <- a[["Total num conserved sites"]]
  }else{
    if (species == "human"){
      #Downlaod Targetscan Species specific site file
      download.file(
        "http://www.targetscan.org/vert_72/vert_72_data_download/Summary_Counts.all_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }else if(any(species %in% c("mouse","rat"))){
      #Downlaod Targetscan Species specific site file
      download.file(
        "http://www.targetscan.org/mmu_72/mmu_72_data_download/Summary_Counts.all_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }
    a$sites <- a[["Total num conserved sites"]]+a[["Total num nonconserved sites"]]
  }
  unlink(tmp)
  a <- as.data.frame(a[a$`Species ID` == spec,])
  a$`Transcript ID` <- gsub("\\..*","",a$`Transcript ID`)
  keep <- c("miRNA family", "Gene Symbol", "sites", "Cumulative weighted context++ score")
  tmp <- a[!duplicated(a[,1:2]),1:2]
  syn <- as.character(tmp[,2])
  names(syn) <- tmp[,1]
  syn <- c(syn, .ens2symbol(species))
  if(keepMers) keep <- c(keep, grep("^Number of", colnames(a), value=TRUE))
  a <- a[,keep]
  colnames(a)[1:4] <- c("set","feature","sites","score")
  a[,1] <- as.factor(a[,1])
  a[,2] <- as.factor(a[,2])
  a[,3] <- as.integer(a[,3])
  a[[4]][a[[4]] == "NULL"] <- 0
  a[,4] <- as.numeric(a[,4])
  #filter for the miRNA species
  a <- DataFrame(a)
  metadata(a)$feature.synonyms <- syn
  metadata(a)$families <- fams[["Seed.m8"]]
  names(metadata(a)$families) <- fams[["MiRBase.ID"]]
  metadata(a)$families <- droplevels(metadata(a)$families)
  a
}


#' importFrom Matrix sparseMatrix miRTarBase
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

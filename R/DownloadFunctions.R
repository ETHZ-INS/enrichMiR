
#' Downlaod Targetscan Position Files
.getTargetScanPos <- function(species = c("human","mouse","rat"), incl_nonconsites = TRUE) {
  species <- match.arg(species)
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116, 
                  stop("No matched species"))
  
  fams <- .fetch_Mirfamilies(species)
  
  # download TargetScan miRNA Positions
  tmp <- tempfile()

  # download TargetScan all conserved miRNA sites positions
  if (species == "human"){
    download.file(
      "http://www.targetscan.org/vert_72/vert_72_data_download/Conserved_Site_Context_Scores.txt.zip", tmp)
    a <- fread(unzip(file.path(tmp)),drop = c("context++ score percentile","weighted context++ score percentile"))
  }else if(any(species %in% c("mouse","rat"))){
    download.file(
      "http://www.targetscan.org/mmu_72/mmu_72_data_download/Conserved_Site_Context_Scores.txt.zip", tmp)
    a <- fread(unzip(file.path(tmp)),drop = c("context++ score percentile","weighted context++ score percentile"))
  }
  a <- merge(a,fams[,c("Seed.m8","MiRBase.ID")], by.x = "miRNA",by.y = "MiRBase.ID")
  a <- a[,.(score=min(`context++ score`), weighted_score = min(`weighted context++ score`)), by = c("Transcript ID","Gene Symbol","Site Type","UTR_start","UTR end","Seed.m8")]
  unlink(tmp)
  
  # if specified download as well the non conserved sites and rbind
  if(incl_nonconsites){
    tmp <- tempfile()
    if (species == "human"){
      download.file(
        "http://www.targetscan.org/vert_72/vert_72_data_download/Nonconserved_Site_Context_Scores.txt.zip", tmp)
      b <- fread(unzip(file.path(tmp)),drop = c("context++ score percentile","weighted context++ score percentile"))
    }else if(any(species %in% c("mouse","rat"))){
      download.file(
        "http://www.targetscan.org/mmu_72/mmu_72_data_download/Nonconserved_Site_Context_Scores.txt.zip", tmp)
      b <- fread(unzip(file.path(tmp)),drop = c("context++ score percentile","weighted context++ score percentile"))
    }
    b <- merge(b,fams[,c("Seed.m8","MiRBase.ID")], by.x = "miRNA",by.y = "MiRBase.ID")
    b <- b[,.(score=min(`context++ score`), weighted_score = min(`weighted context++ score`)), by = c("Transcript ID","Gene Symbol","Site Type","UTR_start","UTR end","Seed.m8")]
    a <- rbind(a,b)
    unlink(tmp) 
  }
  a$`Transcript ID` <- gsub("\\..*","",a$`Transcript ID`)
  a <- as.data.frame(a)
  tmp <- a[!duplicated(a[,1:2]),1:2]
  syn <- as.character(tmp[,2])
  names(syn) <- tmp[,1]
  syn <- c(syn, .ens2symbol(species))
  keep <- c("Seed.m8", "Gene Symbol","UTR_start","UTR end","score","weighted_score")
  a <- a[,keep]
  colnames(a)[1:6] <- c("set","feature","start","end","score","weighted_score")
  a[,1] <- as.factor(a[,1])
  a[,2] <- as.factor(a[,2])
  a[,3] <- as.integer(a[,3])
  a[,4] <- as.integer(a[,4])
  a[[5]][a[[5]] == "NULL"] <- 0
  a[,5] <- as.numeric(a[,5])
  a[[6]][a[[6]] == "NULL"] <- 0
  a[,6] <- as.numeric(a[,6])
  
  #attach the metadata
  a <- DataFrame(a)
  metadata(a)$feature.synonyms <- syn
  metadata(a)$families <- fams[["Seed.m8"]]
  names(metadata(a)$families) <- fams[["MiRBase.ID"]]
  metadata(a)$families <- droplevels(metadata(a)$families)
  a
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
  
  
  tmp <- tempfile()
  if(type=="conserved"){
    # download TargetScan conserved miRNA sites
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
    #Downlaod Targetscan Species specific site file (All Sites)
    if (species == "human"){
      download.file(
        "http://www.targetscan.org/vert_72/vert_72_data_download/Summary_Counts.all_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }else if(any(species %in% c("mouse","rat"))){
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

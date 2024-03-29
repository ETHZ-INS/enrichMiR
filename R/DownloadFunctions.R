# Download Targetscan8 miRNA families
.getTargetscan_miRfamilies <- function(species = c("human","mouse","rat","fish","worm","fly")) {
  species <- match.arg(species)
  
  if(species %in% c("human","mouse","rat")){
    # assign species ID
    spec <- switch( species,
                    human = 9606,
                    mouse = 10090,
                    rat = 10116, 
                    stop("No matched species"))
    
    # download TargetScan miRNA families
    tmp <- tempfile()
    download.file("http://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip", tmp)
    miRFam <- read.delim(unzip(tmp), header=TRUE)
    unlink(tmp)
    
    #filter for the miRNA species
    miRFam <- miRFam[miRFam$Species.ID == spec,]
  }else{
    tmp <- tempfile()
    switch( species,
            fish = download.file("http://www.targetscan.org//fish_62//fish_62_data_download/miR_Family_Info.txt.zip", tmp),
            worm = download.file("http://www.targetscan.org/worm_52/worm_52_data_download/miR_Family_Info.txt.zip", tmp),
            fly = download.file("http://www.targetscan.org/fly_72/fly_72_data_download/miR_Family_Info.txt.zip", tmp),
            stop("No matched species"))
    miRFam <- read.delim(unzip(tmp), header=TRUE)
    unlink(tmp)
  }
  miRFam
}





# Downlaod Targetscan8 Site Position Files
.getTargetScanSites <- function(species = c("human","mouse","rat"), incl_nonconsites = TRUE) {
  species <- match.arg(species)
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116, 
                  stop("No matched species"))
  
  fams <- .getTargetscan_miRfamilies(species)
  
  # download TargetScan miRNA Positions
  tmp <- tempfile()

  # download TargetScan all conserved miRNA sites positions
  if (species == "human"){
    download.file(
      "http://www.targetscan.org/vert_80/vert_80_data_download/Conserved_Site_Context_Scores.txt.zip", tmp)
    a <- fread(unzip(file.path(tmp)),drop = c("context++ score percentile","weighted context++ score percentile"))
  }else if(any(species %in% c("mouse","rat"))){
    download.file(
      "http://www.targetscan.org/mmu_80/mmu_80_data_download/Conserved_Site_Context_Scores.txt.zip", tmp)
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
        "http://www.targetscan.org/vert_80/vert_80_data_download/Nonconserved_Site_Context_Scores.txt.zip", tmp)
      b <- fread(unzip(file.path(tmp)),drop = c("context++ score percentile","weighted context++ score percentile"))
    }else if(any(species %in% c("mouse","rat"))){
      download.file(
        "http://www.targetscan.org/mmu_80/mmu_80_data_download/Nonconserved_Site_Context_Scores.txt.zip", tmp)
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
  
  # add SiteType information
  sit_ty <- data.frame("Site Type" = c(1,2,3,4),"type" = c("7mer-1a","7mer-m8","8mer","6mer"),check.names = FALSE)
  a <- merge(a,sit_ty,by="Site Type")
  
  syn <- c(syn, .ens2symbol(species,level = "gene",report.element = "sym"))
  keep <- c("Seed.m8", "Gene Symbol","UTR_start","UTR end","score","weighted_score","type")
  a <- a[,keep]
  colnames(a)[1:7] <- c("set","feature","start","end","score","weighted_score","type")
  a[,1] <- as.factor(a[,1])
  a[,2] <- as.factor(a[,2])
  a[,3] <- as.integer(a[,3])
  a[,4] <- as.integer(a[,4])
  a[[5]][a[[5]] == "NULL"] <- 0
  a[,5] <- as.numeric(a[,5])
  a[[6]][a[[6]] == "NULL"] <- 0
  a[,6] <- as.numeric(a[,6])
  a[,7] <- as.factor(a[,7])
  
  #attach the metadata
  a <- DataFrame(a)
  metadata(a)$feature.synonyms <- syn
  metadata(a)$families <- fams[["Seed.m8"]]
  names(metadata(a)$families) <- fams[["MiRBase.ID"]]
  metadata(a)$families <- droplevels(as.factor(metadata(a)$families))
  a
  }



# Downlaod Targetscan8 Transcript / Gene Predictions
.getTargetScanPred <- function(species = c("human","mouse","rat","fly"), type=c("conserved","all"), keepMers=FALSE){
  type <- match.arg(type)
  species <- match.arg(species)
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116, 
                  fly = 7227,
                  stop("No matched species"))
  
 
  fams <- .getTargetscan_miRfamilies(species)

  
  tmp <- tempfile()
  if(type=="conserved"){
    # download TargetScan conserved miRNA sites
    if (species == "human"){
      #Downlaod Targetscan file
      download.file(
        "http://www.targetscan.org/vert_80/vert_80_data_download/Summary_Counts.default_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }else if(any(species %in% c("mouse","rat"))){
      #Downlaod Targetscan file
      download.file(
        "http://www.targetscan.org/mmu_80/mmu_80_data_download/Summary_Counts.default_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }else if(species == "fly"){
      #Downlaod Targetscan file
      download.file(
        "http://www.targetscan.org/fly_72/fly_72_data_download/Summary_Counts.default_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }
    a$sites <- a[["Total num conserved sites"]]
  }else{
    #Downlaod Targetscan Species specific file (All Sites)
    if (species == "human"){
      download.file(
        "http://www.targetscan.org/vert_80/vert_80_data_download/Summary_Counts.all_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }else if(any(species %in% c("mouse","rat"))){
      download.file(
        "http://www.targetscan.org/mmu_80/mmu_80_data_download/Summary_Counts.all_predictions.txt.zip", tmp)
      a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    }else if(species == "fly"){
      download.file(
        "http://www.targetscan.org/fly_72/fly_72_data_download/Summary_Counts.all_predictions.txt.zip", tmp)
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
  syn <- c(syn, .ens2symbol(species,level = "gene",report.element = "sym"))
  if(keepMers){
    keep <- c(keep, grep("^Number of", colnames(a), value=TRUE))
    if(type=="conserved") keep <- keep[-grep("nonconserved",keep)]
    if(species != "fly"){
      keep <- keep[-grep("6mer",keep)]
    }}
  a <- a[,keep]
  colnames(a)[1:4] <- c("set","feature","sites","score")
  a[,1] <- as.factor(a[,1])
  a[,2] <- as.factor(a[,2])
  a[,3] <- as.integer(a[,3])
  a[[4]][a[[4]] == "NULL"] <- 0
  a[,4] <- as.numeric(a[,4])
  if(keepMers){
    if(type=="conserved"){
      colnames(a)[grep("8mer",colnames(a))] <- "Sites_8mer"
      colnames(a)[grep("7mer-m8",colnames(a))] <- "Sites_7mer_m8"
      colnames(a)[grep("7mer-1a",colnames(a))] <- "Sites_7mer_1a"
    }else{
      a$`Sites_8mer` <- rowSums(a[,grep("8mer",colnames(a))],na.rm = FALSE)
      a$`Sites_7mer_m8` <- rowSums(a[,grep("7mer-m8",colnames(a))],na.rm = FALSE)
      a$`Sites_7mer_1a` <- rowSums(a[,grep("7mer-1a",colnames(a))],na.rm = FALSE)
      a <- a[,-grep("Number",colnames(a))]
    }
  }
  a <- DataFrame(a)
  metadata(a)$feature.synonyms <- syn
  metadata(a)$families <- fams[["Seed.m8"]]
  names(metadata(a)$families) <- fams[["MiRBase.ID"]]
  metadata(a)$families <- droplevels(as.factor(metadata(a)$families))
  a
}


# Downlaod Targetscan8 Transcript / Gene Predictions for fish and worm
.getTargetScanPred2 <- function(species = c("fish","worm"), type=c("conserved","all"), keepMers=FALSE){
  type <- match.arg(type)
  species <- match.arg(species)
  
  if(species == "fish") message("There is no distinction of 'conserved' & 'all' binding sites in
                                the fish TargetScan database.")
  
  # assign species ID
  spec <- switch( species,
                  fish = 7955,
                  worm = 6239,
                  stop("No matched species"))
  
  
  fams <- .getTargetscan_miRfamilies(species)
  
  
  tmp <- tempfile()
  if(species=="fish"){
    download.file(
        "http://www.targetscan.org//fish_62//fish_62_data_download/Summary_Counts.txt.zip", tmp)
    a <- fread(unzip(file.path(tmp)))
    a$sites <- a[["Total number of sites"]]
    a$`Transcript ID` <- gsub("\\..*","",a$`Transcript ID`)
  }else if(species == "worm"){
    download.file(
        "http://www.targetscan.org/worm_52/worm_52_data_download/Summary_Counts.txt.zip", tmp)
    a <- fread(unzip(file.path(tmp)),drop = c("Aggregate PCT"))
    a$`UTR ID` <- gsub("\\..*","",a$`UTR ID`)
    
    #differentiate between conserved and non-conserved sites
    if(type == "conserved"){
      a <- a[a$`Total num conserved sites` > 0,]
      a$sites <- a$`Total num conserved sites`
    }else if(type == "all"){
      a$sites <- a$`Total num conserved sites` + a$`Total num nonconserved sites`
    }}
  unlink(tmp)
  a <- as.data.frame(a[a$`Species ID` == spec,])
  keep <- c("miRNA family", "Gene Symbol", "sites", "Total context score")
  tmp <- a[!duplicated(a[,1:2]),1:2]
  syn <- as.character(tmp[,2])
  names(syn) <- tmp[,1]
  syn <- c(syn, .ens2symbol(species,level = "gene",report.element = "sym"))
  if(keepMers){
    keep <- c(keep, grep("^Number of", colnames(a), value=TRUE))
    if(type=="conserved" && species == "worm") keep <- keep[-grep("nonconserved",keep)]
    }
  a <- a[,keep]
  colnames(a)[1:4] <- c("set","feature","sites","score")
  a[,1] <- as.factor(a[,1])
  a[,2] <- as.factor(a[,2])
  a[,3] <- as.integer(a[,3])
  a[[4]][a[[4]] == "NULL"] <- 0
  a[,4] <- as.numeric(a[,4])
  if(keepMers){
    if(species == "fish" || type=="conserved"){
      colnames(a)[grep("8mer",colnames(a))] <- "Sites_8mer"
      colnames(a)[grep("7mer-m8",colnames(a))] <- "Sites_7mer_m8"
      colnames(a)[grep("7mer-1a",colnames(a))] <- "Sites_7mer_1a"
    }else{
      a$`Sites_8mer` <- rowSums(a[,grep("8mer",colnames(a))],na.rm = FALSE)
      a$`Sites_7mer_m8` <- rowSums(a[,grep("7mer-m8",colnames(a))],na.rm = FALSE)
      a$`Sites_7mer_1a` <- rowSums(a[,grep("7mer-1a",colnames(a))],na.rm = FALSE)
      a <- a[,-grep("Number",colnames(a))]
    }
  }
  a <- DataFrame(a)
  metadata(a)$feature.synonyms <- syn
  metadata(a)$families <- fams[["Seed.m8"]]
  names(metadata(a)$families) <- fams[["MiRBase.ID"]]
  metadata(a)$families <- droplevels(as.factor(metadata(a)$families))
  a
}


# importFrom Matrix sparseMatrix miRTarBase
.fetch_mirtarbase <- function(species, returnType=c("DataFrame","matrix")){
  options(timeout=500)
  returnType <- match.arg(returnType)
  spec <- switch(species, human="hsa", mouse="mmu", rat="rno", species)
  tmp <- tempfile()
  download.file(paste0("https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/8.0/", 
                       spec, "_MTI.xls", ifelse(spec=="hsa", "x","")),
                destfile=tmp, extra=c(timeout=999))
  e <- readxl::read_excel(tmp)
  e <- e[,c(2,4)]
  e <- as.data.frame(e)
  e[,1] <- as.factor(e[,1])
  e[,2] <- as.factor(e[,2])
  e <- e[!duplicated(e),]
  if(returnType=="DataFrame"){
    e <- DataFrame(e)
    colnames(e) <- c("set","feature")
    syn <- .ens2symbol(species,level = "gene",report.element = "sym")
    metadata(e)$feature.synonyms <- syn
    return(e)
  }
  sparseMatrix( as.numeric(e[,2]), as.numeric(e[,1]), 
                dimnames=list(levels(e[,2]), levels(e[,1])) )
}




#Downlaod RBP position info from the ornament database
.getOrnamentDB <- function(species = c("human","mouse"), type=c("coding","non-coding"),MSSthreshold = "50%"){
  message("This function will take quite some time")
  species <- match.arg(species)
  type <- match.arg(type)

  #Download raw data by species and type
  tmp <- tempfile()
  options(timeout=750)
  if(species == "human"){
    if(type == "coding"){
      download.file(
        "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Homo_sapiens_cDNA_oRNAment.csv.gz", tmp)
    }else if(type == "non-coding"){
      download.file(
        "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Homo_sapiens_ncRNA_oRNAment.csv.gz", tmp)
    }else{stop("Choose either 'coding' or 'non-coding'")}
  }else if(species == "mouse"){
    if(type == "coding"){
      download.file(
        "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Mus_musculus_cDNA_oRNAment.csv.gz", tmp)
    }else if(type == "non-coding"){
      download.file(
        "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Mus_musculus_ncrna_oRNAment.csv.gz", tmp)
    }else{stop("Choose either 'coding' or 'non-coding'")}
  }else{stop("Choose either 'mouse' or 'human'")}
  options(timeout=60)
  #Unzip and read file
  gunzip(file.path(tmp),paste0(file.path(tmp),".csv"), remove=FALSE)
  a <- fread(paste0(file.path(tmp),".csv"),drop = c(1,9,11:12))
  colnames(a) <- c("ensembl_transcript_id_INT", "gene_biotype", "transcript_biotype", "transcript_position", "rbp_code", "score", "unpaired_probability","region")
  unlink(tmp)
  
  # Get Ensembl gene & transcript names
  tmp <- tempfile()
  if(species == "human"){
    download.file(
      "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Homo_sapiens_string_to_int_ID_conversion.csv.gz", tmp)
  }else{
    download.file(
      "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Mus_musculus_string_to_int_ID_conversion.csv.gz", tmp)
  }
  gunzip(file.path(tmp),paste0(file.path(tmp),".csv"), remove=FALSE)
  ens <- fread(paste0(file.path(tmp),".csv"), drop = 5)
  colnames(ens) <- c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "ensembl_transcript_id_INT")
  if(species == "mouse") a$ensembl_transcript_id_INT <- trunc(a$ensembl_transcript_id_INT)
  a <- merge(a,ens,by = "ensembl_transcript_id_INT", all.x = TRUE)
  a$ensembl_transcript_id_INT <- NULL
  unlink(tmp)
  
  
  # Get the tx lenght
  # They use GRCx38 and Ensembl v97
  library(AnnotationHub)
  ah <- AnnotationHub()
  if(species == "human"){
    en <- query(ah, c("EnsDb", "Homo sapiens", "v97"))
  }else{
    en <- query(ah, c("EnsDb", "Mus musculus", "v97"))
  }
  EnsDB <- ah[[en$ah_id]]
  len <- transcriptLengths(EnsDB,with.utr3_len = TRUE)
  len <- as.data.table(len)
  a <- merge(a,len[,c("tx_id","utr3_len","tx_len")],by.x = "ensembl_transcript_id",by.y = "tx_id",all.x = TRUE)
  
  # Filter gene biotype, transcript biotype and get 3'UTR start
  if(type == "coding"){
    #protein coding genes
    a <- a[a$gene_biotype == 1,]
    a <- a[a$transcript_biotype == 1,]
  }else{
    #lincRNA and associated genes
    a <- a[a$gene_biotype == 13,]
    a <- a[a$transcript_biotype == 5,]
    a$region <- NULL
    a$utr3_len <- NULL
  }
  a$gene_biotype <- NULL
  a$transcript_biotype <- NULL
  
  
  
  #Filter for the specified MSS score
  if(MSSthreshold != "50%"){
    MSS <- switch(MSSthreshold,
                  '45%' = "0_45",
                  '40%' = "0_40",
                  '35%' = "0_35",
                  '30%' = "0_30",
                  '25%' = "0_25",
                  '20%' = "0_20",
                  '15%' = "0_15",
                  '10%' = "0_10",
                  '5%' = "0_05",
                  '1%' = "0_01",
                  stop("Possible thresshold values are: '50%','45%','40%','35%','30%','25%','20%','15%','10%','5%','1%'"))
    tmp <- tempfile()
    tmp2 <- tempdir()
    download.file(
      "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/PMM_thresholds.tar.gz", tmp)
    untar(file.path(tmp),exdir = tmp2)
    mss_file <- fread(paste0(tmp2,"/PMM_thresholds/MSS-prime-",MSS,".csv"))
    mss <- mss_file$score
    names(mss) <- mss_file$motif
    idx <- a$score
    names(idx) <- a$rbp_code
    w <- ifelse(idx < mss[names(idx)],FALSE,TRUE)
    a <- a[w,]
    unlink(tmp)
    unlink(tmp2)
  }
  
  # Get the RBP names
  tmp <- tempfile()
  download.file(
    "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/RBP_id_encoding.csv.gz", tmp)
  gunzip(file.path(tmp),paste0(file.path(tmp),".csv"), remove=FALSE)
  rbp <- fread(paste0(file.path(tmp),".csv"))
  colnames(rbp) <- c("rbp_code","RBP")
  a <- merge(a,rbp,by = "rbp_code",all.x = TRUE)
  a$rbp_code <- NULL
  
  # decrease file size
  a$ensembl_transcript_id <- as.factor(a$ensembl_transcript_id) 
  a$score <- as.integer(round(a$score,3)*1000)
  a$unpaired_probability <- as.integer(round(a$unpaired_probability,3)*1000)
  if(!is.null(a$region)) a$region <- as.factor(a$region)
  a$ensembl_gene_id <- as.factor(a$ensembl_gene_id)
  a$external_gene_name <- as.factor(a$external_gene_name)
  if(!is.null(a$utr3_len)) a$utr3_len <- as.integer(a$utr3_len)
  a$tx_len <- as.integer(a$tx_len)
  a$RBP_Motif <- as.factor(a$RBP)
  
  #save as GRanges ?
  # a <- GRanges(seqnames = a$ensembl_transcript_id,ranges = IRanges(start = a$transcript_position,width = 7L), "score" = a$score, "unpaired_probability" = a$unpaired_probability,"region" = a$region, 
  #              "ensembl_gene_id" = a$ensembl_gene_id, "gene_name" = a$external_gene_name, "utr3_len" = a$utr3_len, "tx_len" = a$tx_len, "RBP" = a$RBP)
  a
}



#Import Function
# if tx=FALSE, the longest transcript of each gene gets imported
.import_oRNAment <- function(x, tx = TRUE, only3utr = TRUE, aggregatesites = TRUE){
 
  if(is.null(x$region) && only3utr) stop("The 3'-UTR region can't be exported for lncRNAs")
  
  # filter 3'utr
  if(only3utr){
      x <- x[x$region == "3;3"]
      #get the position relative within in the 3'UTR
      x$transcript_position <- x$transcript_position - (x$tx_len - x$utr3_len)
    }
  
  # only transcript info
  if(isFALSE(tx)){
    if(is.null(x$utr3_len)){
      x <- x[, .SD[tx_len == max(tx_len)], by = .(ensembl_gene_id)]
    }else{
      x <- x[, .SD[utr3_len == max(utr3_len)], by = .(ensembl_gene_id)]
    }
    # check in case there is two transcripts with the same length or utr length
    x <- x[, .SD[ensembl_transcript_id == min(as.character(ensembl_transcript_id))], by = .(ensembl_gene_id)]
    syn <- x$external_gene_name
    names(syn) <- x$ensembl_gene_id
  }
  
  #aggregate all motifs per RBP
  x$RBP <- gsub(" \\(.*","",x$RBP_Motif)
  fams <- unique(x[,c("RBP","RBP_Motif")])

  #aggregate sites
  if(aggregatesites){
    x <- x[, sites:= .N, by = .(ensembl_transcript_id,RBP)]
    if(tx){
      ll = DataFrame("set" = x$RBP, "feature" = x$ensembl_transcript_id, "sites" = x$sites)
      
    }else{
      ll = DataFrame("set" = x$RBP, "feature" = x$external_gene_name, "sites" = x$sites)
      metadata(ll)$feature.synonyms <- syn
    }

  }else{
    if(tx){
      ll = DataFrame("set" = x$RBP, "feature" = x$ensembl_transcript_id, "start" = x$transcript_position, "end" = x$transcript_position + 6)
    }else{
      ll = DataFrame("set" = x$RBP, "feature" = x$external_gene_name, "start" = x$transcript_position, "end" = x$transcript_position + 6)
      metadata(ll)$feature.synonyms <- syn
    }
  }
  metadata(ll)$families <- fams$RBP
  names(metadata(ll)$families) <- fams$RBP_Motif
  metadata(ll)$families <- droplevels(as.factor(metadata(ll)$families))
  ll$set <- as.factor(ll$set)
  ll$feature <- as.factor(ll$feature)
  return(ll)
}



#Prepares a DataFrame that is compatible with EnrichMiR from scanMiR "runFullScan" aggregations (sc)
.prepScanMir <- function(sc,species = c("human","mouse","rat"),level = c("gene","transcript"),include6mers = FALSE, canonicalUTRSites = TRUE){

  #get the Targetscan miRNA families
  species <- match.arg(species)
  level <- match.arg(level)
  fams <- .getTargetscan_miRfamilies(species)
  
  if(class(sc) == "IndexedFst") sc <- as.data.frame(sc)
  if(canonicalUTRSites) sc[colnames(sc) %in% c("non-canonical","ORF.canonical","ORF.non-canonical")] <- NULL
  if(!include6mers) sc[["6mer"]] <- NULL
  if(canonicalUTRSites){
    sc$sites <- rowSums(sc[,grep("mer",colnames(sc))])
  }else{
    sc$sites <- rowSums(sc[,grep("mer|canonical",colnames(sc))])
  }
  sc$repression <- sc$repression / 1000
  colnames(sc)[colnames(sc) == "repression"] <- "score"
  colnames(sc)[colnames(sc) == "miRNA"] <- "set"
  colnames(sc)[colnames(sc) == "transcript"] <- "feature"
  colnames(sc)[colnames(sc) == "8mer"] <- "Sites_8mer"
  colnames(sc)[colnames(sc) == "7mer"] <- "Sites_7mer"
  if(include6mers) colnames(sc)[colnames(sc) == "6mer"] <- "Sites_6mer"
  if(!canonicalUTRSites){
    colnames(sc)[colnames(sc) == "non-canonical"] <- "Non-canonical"
    colnames(sc)[colnames(sc) == "ORF.canonical"] <- "Sites_ORF.canonical"
    colnames(sc)[colnames(sc) == "ORF.non-canonical"] <- "Sites_ORF.non-canonical"
  }
  sc <- sc[sc$sites > 0,]
  
  #filter out non protein-coding biotypes
  anno <- .ensbiotype_df(species = species,level = "tx")
  anno <- anno[anno$transcript_biotype == "protein_coding",]
  sc <- sc[sc$feature %in% anno$ensembl_transcript_id,]
  
  #potential gene aggregation
  if(level == "gene"){
    anno <- .enstx2gene_df(species = species)
    sc <- merge(sc,anno,by.x = "feature",by.y = "ensembl_transcript_id",all.x = TRUE)
    sc <- sc[!is.na(sc$ensembl_gene_id),]
    sc <- sc[order(sc$score),]
    sc <- sc[!duplicated(sc[c("set","ensembl_gene_id")]),]
    sc$feature <- NULL
    colnames(sc)[colnames(sc) == "ensembl_gene_id"] <- "feature"
  }
  
  if(canonicalUTRSites){
    sc <- sc[,c(grep("feature",colnames(sc)),grep("set",colnames(sc)),grep("sites",colnames(sc)),grep("score",colnames(sc)),grep("mer",colnames(sc)))]
  }else{
    sc <- sc[,c(grep("feature",colnames(sc)),grep("set",colnames(sc)),grep("sites",colnames(sc)),grep("score",colnames(sc)),grep("mer",colnames(sc)), grep("canonical",colnames(sc)))]
  }
  sc[,1] <- as.factor(sc[,1])
  sc[,2] <- as.factor(sc[,2])
  sc[,3] <- as.integer(sc[,3])
  sc[,4] <- as.numeric(sc[,4])
  
  
  #add the ensembl/symbol conversion
  if(level == "gene"){
    syn <- .ens2symbol(species = species,level = "gene",report.element = "ens")
  }else{
    syn <- .ens2symbol(species = species,level = "tx",report.element = "ens")
  }
  
  sc <- DataFrame(sc)
  metadata(sc)$feature.synonyms <- syn
  metadata(sc)$families <- fams[["Seed.m8"]]
  names(metadata(sc)$families) <- fams[["MiRBase.ID"]]
  metadata(sc)$families <- droplevels(as.factor(metadata(sc)$families))
  sc
}


.ens2symbol <- function(species=c("human","mouse","rat","fly","worm","fish"),level = c("tx","gene"),
                        report.element=c("sym","ens")){
  library(biomaRt)
  species <- match.arg(species)
  level <- match.arg(level)
  report <- match.arg(report.element)
  symbol <- switch(species,
                   human="hgnc_symbol",
                   mouse="mgi_symbol",
                   rat="rgd_symbol",
                   fly = "external_gene_name",
                   worm = "external_gene_name",
                   fish = "hgnc_symbol",
                   "uniprot_gn_symbol")
  ensembl <- switch(species,
                    human = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl"),
                    mouse = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl"),
                    rat = useMart("ENSEMBL_MART_ENSEMBL",dataset="rnorvegicus_gene_ensembl"),
                    fly = useMart("ENSEMBL_MART_ENSEMBL",dataset="dmelanogaster_gene_ensembl"),
                    worm = useMart("ENSEMBL_MART_ENSEMBL",dataset="celegans_gene_ensembl"),
                    fish = useMart("ENSEMBL_MART_ENSEMBL",dataset="drerio_gene_ensembl"),
                    stop("No matched species"))
  if(level == "gene"){
    anno <- getBM(attributes=c("ensembl_gene_id",symbol), mart = ensembl)
  }else{
    anno <- getBM(attributes=c("ensembl_transcript_id",symbol), mart = ensembl)
  }
  anno <- anno[anno[,2]!="",]
  anno <- anno[!duplicated(anno[,1]),]
  if(report=="sym"){
    e2s <- anno[,2]
    names(e2s) <- anno[,1]
  }else {
    e2s <- anno[,1]
    names(e2s) <- anno[,2]
  }
  e2s
}

# @importFrom biomaRt useMart getBM
.ensbiotype_df <- function(species=c("human","mouse","rat"),level = c("tx","gene")){
  library(biomaRt)
  species <- match.arg(species)
  level <- match.arg(level)
  ensembl <- switch(species,
                    human = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl"),
                    mouse = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl"),
                    rat = useMart("ENSEMBL_MART_ENSEMBL",dataset="rnorvegicus_gene_ensembl"), 
                    stop("No matched species"))
  if(level == "gene"){
    anno <- getBM(attributes=c("ensembl_gene_id","gene_biotype"), mart = ensembl)
  }else{
    anno <- getBM(attributes=c("ensembl_transcript_id","transcript_biotype"), mart = ensembl)
  }
  anno
}

.enstx2gene_df <- function(species=c("human","mouse","rat")){
  library(biomaRt)
  species <- match.arg(species)
  ensembl <- switch(species,
                    human = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl"),
                    mouse = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl"),
                    rat = useMart("ENSEMBL_MART_ENSEMBL",dataset="rnorvegicus_gene_ensembl"), 
                    stop("No matched species"))
  anno <- getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id"), mart = ensembl)
  anno
}



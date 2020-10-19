#' importTargetScan
#'
#' import a targetScan conserved targets file (conserved site and contextScore).
#'
#' @param ts.filepath Path to the targetScan conserved targets file.
#' @param featureIdentifier The type of feature identifier to use as feature name; either "symbol", "transcript", or "geneid".
#' @param recapitalizeGenes Logical; whether to recapitalize genes to the mouse standard (default FALSE)
#'
#' @return a data.frame.
#'
#' @export
importTargetScan <- function(ts.filepath, featureIdentifier="symbol", recapitalizeGenes=FALSE,species = c("human","mouse","rat")){
    featureIdentifier <- match.arg(featureIdentifier, c("symbol","transcript","geneID"))
    ff <- switch( tolower(featureIdentifier),
                  "symbol"="Gene Symbol",
                  "ensembl"="ensembl_gene_id")
    sp <- switch(species,
                 human = "hsapiens_gene_ensembl",
                 mouse = "mmusculus_gene_ensembl",
                 rat = "mmusculus_gene_ensembl", 
                 stop("No matched species"))
    ts2 <- read.csv(ts.filepath,header=T,stringsAsFactors=FALSE)
    if(any(grepl("Total num nonconserved sites",colnames(ts2)))){
        ts2$sites <- ts2$`Total num conserved sites` + ts2$`Total num nonconserved sites`
    }else{
        ts2$sites <- ts2$`Total num conserved sites`
    }
    ts2 <- as.data.table(ts2)
    ts2 <- ts2[ts2$`Cumulative weighted context++ score` != "NULL",]
    ts2$`Cumulative weighted context++ score` <- as.numeric(ts2$`Cumulative weighted context++ score`)
    if(featureIdentifier == "symbol"){
        TS[[ff]] <- switch(species,
                           human = toupper(TS[[ff]]),
                           mouse = paste0(toupper(substring(TS[[ff]],1,1)), tolower(substring(TS[[ff]],2))),
                           rat = paste0(toupper(substring(TS[[ff]],1,1)), tolower(substring(TS[[ff]],2))), 
                           stop("No matched species"))}
    if(featureIdentifier == "ensembl"){
        ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset=sp)
        anno <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id"),mart = ensembl)
        TS <- merge(TS,anno,by.x = "Transcript ID",by.y = "ensembl_transcript_id")
    }
    ts2 <- ts2[,.(sites = max(abs(sites)),score = max(abs(ts2$`Cumulative.weighted.context...score`))),by = list(set = ts2$miRNA.family,feature=ts2[[ff]])]
    ts2 <- DataFrame(ts2)
    return(ts2)
}


#Pay attention, Ensembl mapping with rat does not make much sense!!

.initializeTargetscan <- function(TS,featureIdentifier="symbol",species = c("human","mouse","rat")){
    species <- match.arg(species)
    featureIdentifier <- match.arg(featureIdentifier, c("symbol","ensembl"))
    ff <- switch( tolower(featureIdentifier),
                  "symbol"="Gene Symbol",
                  "ensembl"="ensembl_gene_id")
    sp <- switch(species,
                       human = "hsapiens_gene_ensembl",
                       mouse = "mmusculus_gene_ensembl",
                       rat = "mmusculus_gene_ensembl", 
                       stop("No matched species"))
    TS <- as.data.table(TS)
    if(any(grepl("Total num nonconserved sites",colnames(TS)))){
    TS$sites <- TS$`Total num conserved sites` + TS$`Total num nonconserved sites`
    }else{
    TS$sites <- TS$`Total num conserved sites`
    }
    TS <- TS[TS$`Cumulative weighted context++ score` != "NULL",]
    TS$`Cumulative weighted context++ score` <- as.numeric(TS$`Cumulative weighted context++ score`)
    if(featureIdentifier == "symbol"){
        TS[[ff]] <- switch(species,
                       human = toupper(TS[[ff]]),
                       mouse = paste0(toupper(substring(TS[[ff]],1,1)), tolower(substring(TS[[ff]],2))),
                       rat = paste0(toupper(substring(TS[[ff]],1,1)), tolower(substring(TS[[ff]],2))), 
                       stop("No matched species"))}
    if(featureIdentifier == "ensembl"){
        ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset=sp)
        anno <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id"),mart = ensembl)
        TS <- merge(TS,anno,by.x = "Transcript ID",by.y = "ensembl_transcript_id")
    }
    TS <- TS[,.(sites = max(abs(sites)),score = max(abs(`Cumulative weighted context++ score`))),by = list(set = TS$`miRNA family`,feature=TS[[ff]])]
    TS <- DataFrame(TS)
    TS
}



.plMergeList <- function(ll,...){
  if(!is.null(names(ll))) for(i in 1:length(ll)) colnames(ll[[i]]) <- paste(names(ll)[i],colnames(ll[[i]]),sep=".")
  while(length(ll)>1){
    x <- merge(ll[[1]],ll[[2]],by="row.names",...)
    row.names(x) <- x[,1]
    ll[[1]] <- x[,-1]
    ll[[2]] <- NULL
  }
  ll[[1]]
}

#' recapitalizeGenes
#'
#' A utility function to reformat gene names.
#'
#' @param x A character vector or an array (in which case the function will be applied to row.names)
#' @param gformat The gene format, either 'human' (all caps) or 'mouse' (first letter capitalized)
#' 
#' @return An object of the same type and dimensions as `x`.
#'
#' @export
recapitalizeGenes <- function(x, gformat="mouse"){
  gformat <- match.arg(gformat, choices=c("human","mouse"))
  if(is(x,'data.frame') | is(x,'matrix')){
    row.names(x) <- recapitalizeGenes(row.names(x), gformat)
    return(x)
  }
  switch(gformat,
    mouse=sapply(x,FUN=function(x){ paste0(toupper(substring(x,1,1)), tolower(substring(x,2))) }),
    human=sapply(x,toupper)
  )
}

#' recapitalizeMiRs
#'
#' A utility function to reformat miRNA names so that there are no capital letters except the 'miR-'.
#'
#' @param x A character vector or an array (in which case the function will be applied to row.names)
#' 
#' @return An object of the same type and dimensions as `x`.
#'
#' @export
recapitalizeMiRs <- function(x){
  if(is(x,'data.frame') | is(x,'matrix')){
    if(is.null(row.names(x))) stop("Given an array without row.names... please apply to the column containing the miRNA names.")
    row.names(x) <- recapitalizeMiRs(row.names(x))
    return(x)
  }
  x <- gsub("mir-","miR-",tolower(x))
}

.cleanMiRname <- function(x){ paste(strsplit(x,"-",fixed=T)[[1]][-1],collapse="-") }


.dea2binary <- function( dea, th=0.05, th.alfc=0, min.at.th=20, alt.top=50, 
                         restrictSign=NULL, verbose=TRUE ){
  dea <- .homogenizeDEA(dea)
  if(!is.null(restrictSign)){
    if(!(restrictSign %in% c(-1,1)))
      stop("`restrictSign`, if given, should be -1 or 1.")
    dea$FDR[sign(dea$logFC)!=restrictSign] <- 1
  }
  dea <- dea[order(dea$FDR, dea$PValue),]
  if(sum(bi <- (dea$FDR<=th & abs(dea$logFC)>=th.alfc),na.rm=TRUE) < min.at.th){
    if(verbose) message("Insufficient genes passing the defined FDR; will use",
                        "the top ", alt.top, " genes.")
    x <- rep(c(TRUE,FALSE), c(alt.top,nrow(dea)-alt.top))
  }else{
    x <- bi
  }
  names(x) <- row.names(dea)
  x
}


.dea2sig <- function( dea, field=NULL ){
  if(is.null(field)){
    dea <- .homogenizeDEA(dea)
    x <- sign(dea$logFC)*-log10(dea$FDR)
  }else{
    x <- dea[[field]]
  }
  names(x) <- row.names(dea)
  x
}

.homogenizeDEA <- function(x, keepTop=TRUE){
  if(is(x,"data.table")){
    if(any(duplicated(x[[1]]))){
      if(keepTop){
        x <- x[order(x[[head(grep("padj|adj\\.P\\.Val|q_value|qval", colnames(x)),1)]]),]
        x <- x[!duplicated(x[[1]]),]
      }else{
        x <- aggregate(x[,-1,drop=FALSE], by=list(gene=x[[1]]), FUN=mean)
      }
    }
    x <- data.frame(row.names=x[[1]], as.data.frame(x[,-1,drop=FALSE]))
  }
  x <- as.data.frame(x)
  w <- grep("^ENS",row.names(x))
  row.names(x)[w] <- gsub("\\..*","",row.names(x)[w])

  colnames(x) <- gsub("log2FoldChange|log2Fold|log2FC|log2\\(fold_change\\)|log2\\.fold_change\\.",
                      "logFC", colnames(x))
  
  abf <- colnames(df)[which(colnames(df) %in% c("meanExpr", "AveExpr", 
                                                "baseMean", "logCPM"))]
  if (length(abf) == 1) {
    x$meanExpr <- df[, abf]
    if (abf == "baseMean") 
      x$meanExpr <- log(x$meanExpr + 1)
  }else if(all(c("value_1","value_2") %in% colnames(x))){ # cufflinks
    x$meanExpr <- log(1+x$value_1+x$value_2)
  }
  colnames(x) <- gsub("P\\.Value|pvalue|p_value|pval", "PValue", colnames(x))
  colnames(x) <- gsub("padj|adj\\.P\\.Val|q_value|qval", "FDR", colnames(x))
  if (!("FDR" %in% colnames(x))) 
    x$FDR <- p.adjust(x$PValue, method = "fdr")
  f <- grep("^logFC$",colnames(x),value=TRUE)
  if(length(f)==0) f <- grep("logFC",colnames(x),value=TRUE)
  if(length(f)==0) warning("No logFC found.")
  if(length(f)>1){
    message("Using ",f[1])
    x[["logFC"]] <- x[[f[1]]]
  }
  x
}


#' dround
#'
#' Trim to a certain number of digits (similar to `format(...,digits=digits)`, 
#' except that the output remains numeric)
#'
#' @param x A vector of numeric values, or a data.frame or similar
#' @param digits The number of digits to keep
#' @param roundGreaterThan1 Whether to trim also numbers greater than 1 (default FALSE)
#'
#' @return A object of same dimensions as `x`
#' @export
#'
#' @examples
#' dround( c(0.00002345, 554356, 12.56) )
dround <- function(x, digits=3, roundGreaterThan1=FALSE){
  if(!is.null(dim(x))){
    for(i in seq_len(ncol(x))){
      if(!is.integer(x[,i]) && is.numeric(x[,i])){
        tryCatch({
          x[,i] <- dround(as.numeric(x[,i]), digits, roundGreaterThan1)
        }, error=function(e) warning(e))
      }
    }
    return(x)
  }  
  if(roundGreaterThan1){
    w <- 1:length(x)
  }else{
    w <- which(abs(x)<1)
  }
  if(length(w)==0) return(x)
  e <- ceiling(-log10(abs(x[w])))
  y <- x
  y[w] <- round(10^e*x[w],digits-1)/10^e
  y[x==0] <- 0
  y
}



.agDF <- function(df, new.rn, 
                  match.names=(!is.null(names(new.rn)) && length(df)!=length(new.rn))){
  if(match.names){
    names(a) <- a <- row.names(df)
    a[names(new.rn)] <- as.character(new.rn)
    new.rn <- a[row.names(df)]
  }
  tt <- split(names(new.rn),new.rn)
  tt <- data.frame(row.names=names(tt),
                   members=sapply(tt, FUN=function(x) paste(sort(x),collapse=";")))
  if(ncol(df)==0) return(data.frame(row.names=unique(new.rn)))
  if(!any(duplicated(new.rn))){
    row.names(df) <- new.rn
    df <- merge(df,tt,by = 0, all.x = TRUE)
    row.names(df) <- df$Row.names
    df[,-c(1),drop=FALSE]
  }
  df <- aggregate(df, by=list(RN=new.rn), FUN=function(x){
    if(is.factor(x)) x <- as.character(x)
    if(length(x)==1) return(x)
    if(is.numeric(x)) return(x[which.max(abs(x))])
    paste(sort(x), collapse=";")
  })
  row.names(df) <- df[,1]
  df <- merge(df,tt,by = 0, all.x = TRUE)
  row.names(df) <- df$Row.names
  df[,-c(1,2),drop=FALSE]
}

.filterTranscripts <- function(x, minProp=0.9, minLogCPM=1){
  if(!all(c("transcript","gene","logCPM") %in% colnames(x))) stop("Malformed input.")
  gs <- rowsum(exp(x$logCPM), x$gene)
  x[ (exp(x$logCPM)/gs[x$gene,1]) > minProp & x$logCPM>minLogCPM, ]
}

# triggers an error if the sets are not formatted correctly, and return a 
# logical indicating whether the sets are in matrix format
.checkSets <- function(sets, requiredColumns=c(), matrixAlternative=FALSE){
  if( !isFALSE(matrixAlternative) && 
      (is.matrix(sets) || is(sets,"sparseMatrix")) ){
    matrixAlternative <- match.arg(matrixAlternative, c("logical","numeric"))
    if(is(sets,"sparseMatrix")){
      if(matrixAlternative=="logical") stopifnot(is(sets,"lgCMatrix"))
      if(matrixAlternative=="numeric") stopifnot(is(sets,"dgCMatrix"))
    }else{
      if(matrixAlternative=="logical") stopifnot(is.logical(sets))
      if(matrixAlternative=="numeric") stopifnot(is.numeric(sets))
    }
    return(TRUE)
  }
  stopifnot(is.data.frame(sets) || is(sets, "DFrame"))
  stopifnot(c("feature","set",requiredColumns) %in% colnames(sets))
  FALSE
}

.is.matrix <- function(x){
  is.matrix(x) || is(x,"sparseMatrix") || is(x,"DelayedArray")
}

#' @import S4Vectors
.list2DF <- function(sets){
  if(!is.null(dim(sets)) && !all(c("feature","set") %in% colnames(sets)))
    stop("Malformed `sets`.")
  if(is(sets, "DataFrame")) return(sets)
  if(is.data.frame(sets)){
    at <- attributes(sets)
    sets <- DataFrame(sets)
    for(f in intersect(at, c("sets.properties","feature.synonyms")))
      metadata(sets)[[f]] <- at[[f]]
    return(sets)
  }
  if(is.null(names(sets))) stop("The sets should be named!")
  if(is.list(sets[[1]])){
    if(all(names(sets[[1]]) %in% c("tfmode","likelihood"))){
      # regulon object
      sets <- lapply(sets, FUN=function(x){
        data.frame(feature=names(x[[1]]), score=x$tfmode*x$likelihood)
      })
      return(cbind(set=rep(names(sets),sapply(sets,nrow)),
                   do.call(rbind, sets)))
    }else{
      stop("`sets` has an unknown format")
    }
  }
  y <- DataFrame( set=factor(rep(sets, lengths(sets))) )
  if(!is.null(names(sets[[1]])) && is.numeric(sets[[1]])){
    y$feature <- unlist(lapply(sets, names))
    y$score <- unlist(lapply(sets, as.numeric))
  }else{
    y$feature <- unlist(lapply(sets, as.character))
  }
  y$feature <- as.factor(y$feature)
  y
}




#' getHumanMirExp
#'
#' Get the human miRNA expression in a given tissue or celltype from the 
#' `microRNAome` package 
#' (\href{http://genome.cshlp.org/content/27/10/1769}{McCall et al., 2017}).
#'
#' @param x Desired tissue or celltype
#'
#' @return If `x` is NULL, returns a vector of the possible values. If `x` is 
#' given and matches a tissue/celltype of the dataset, returns the average miRNA
#' expression as logCPM, i.e. log(1+counts per million).
#' @export
#'
#' @examples
#' head(getHumanMirExp("thyroid"))
getHumanMirExp <- function(x=NULL){
  data("microRNAome", package="microRNAome")
  if(is.null(x)) return(sort(unique(microRNAome$cell_tissue)))
  if(!(x %in% microRNAome$cell_tissue)) return(NULL)
  x <- assay(microRNAome)[,microRNAome$cell_tissue==x,drop=FALSE]
  row.names(x) <- gsub("/.+","",(row.names(x)))
  x <- rowSums(x)
  sort(log1p(10^6 * x/sum(x)),decreasing=TRUE)
}

#' getMouseMirExp
#'
#' Get the mouse miRNA expression in a given tissue or celltype, based on the
#' GSE119661 dataset 
#' (\href{https://doi.org/10.1093/nar/gkaa323}{Kern et al., 2020}) 
#' supplemented with some celltypes from the GSE30286 dataset 
#' (\href{https://doi.org/10.1016/j.neuron.2011.11.010}{He et al., 2012}).
#'
#' @param x Desired tissue or celltype
#'
#' @return If `x` is NULL, returns a vector of the possible values. If `x` is 
#' given and matches a tissue/celltype of the dataset, returns the average miRNA
#' expression as logCPM, i.e. log(1+counts per million).
#' @export
#'
#' @examples
#' head(getMouseMirExp("Muscle"))
getMouseMirExp <- function(x=NULL){
  data("miRNAexpMouse", package="enrichMiR")
  if(is.null(x)) return(sort(colnames(miRNAexpMouse)))
  if(!(x %in% colnames(miRNAexpMouse))) return(NULL)
  sort(miRNAexpMouse[,x], decreasing=TRUE)
}


.exampleBackground <- function(){
  c("PLEKHB2", "FBXO21", "PIGS", "SLC7A1", "YKT6", "RABL6", "SLC1A5", 
    "OCRL", "SH3RF1", "SLC9A1", "SLC7A5", "SRPRA", "G6PC3", "ARL2", 
    "CTDSPL", "CSRP1", "TBC1D22B", "ZER1", "SLC52A2", "GNPDA1", "PRDM4", 
    "VAMP3", "TMEM87A", "SEPTIN2", "SHCBP1", "RELL1", "PRTFDC1", 
    "CPEB1", "PHKA1", "RNF38", "TMED8", "CASP7", "PROSER1", "LMNB2", 
    "DSTYK", "TMEM216", "RAVER1", "PRUNE1", "ASF1B", "NCDN", "TGFBRAP1", 
    "CS", "MTHFD2", "SFT2D1", "ARHGAP1", "IQGAP1", "ATN1", "CTDNEP1", 
    "VPS37B", "PLCD3", "PKM", "POLR3D", "SLC25A6", "PRKRA", "HECTD3", 
    "SULF1", "SERINC5", "DYNC1LI2", "BCAT2", "VPS4A", "AL356776.1", 
    "DMAC1", "MROH8", "MSH2", "PANX2", "PAQR9", "SWI5", "DHX9P1", 
    "RAP1GDS1", "OLFM2", "PGF", "TBC1D12", "TNFSF12", "ANP32BP1", 
    "OSGIN2", "LMCD1", "CNOT9", "TTC6", "YTHDF2", "BOD1", "ZMYM3", 
    "USF1", "STK4", "PYCR3", "AC104083.1", "NSUN5P2", "PTAFR", "CFAP58", 
    "CDKN2AIP", "ARPC2", "HMGN2P46", "PJA1", "HSPA6", "AC079416.3", 
    "MYH9", "MRTFB", "SDK1", "HNRNPA1P69", "RAB10", "BUD13", "NDUFA2", 
    "STK31", "SEC11B", "ROMO1", "KDELC1P1", "EXOSC5", "AC092807.3", 
    "AC004492.1", "PRSS53", "MTRNR2L6", "RPS2P6", "RPL5P1", "AC013391.2", 
    "PNPT1P1", "CALM2", "AC027801.2", "ANO8", "RAB2B", "AL136038.4", 
    "PPM1B", "BHLHE41", "SBF1", "LIFR-AS1", "PSMB4", "POGK", "CSDC2", 
    "PPIAP26", "PHRF1", "NPM1P9", "PDIA2", "BCRP2", "GAPDHP58", "CLEC16A", 
    "AC234775.3", "EP400", "VTI1B", "PTPN6", "SIMC1", "AL034397.1", 
    "RBBP8", "PEX14", "ABCA11P", "TEN1-CDK3", "ZYG11B", "AGO1", "TAF3", 
    "GBX2", "PLIN5", "PKDCC", "SMNDC1", "RUSC2", "FAM20C", "HBP1", 
    "SDHAF1", "ZBTB40-IT1", "TMEM14C", "TUFM", "SH3GL2", "MDM4", 
    "MIR3936HG", "FAM210A", "SNRPGP10", "CCAR1", "RORB", "AL645940.1", 
    "RPN1", "LSM6", "PPRC1", "CD70", "AC027307.3", "ZNF789", "BRSK2", 
    "MEGF6", "DPYSL3", "CRYBG3", "DMAP1", "TMEM101", "CDKN1B", "UCP1", 
    "ZNF222", "LINC00926", "RGS7", "EML4", "PCBP1", "BTF3P8", "INTS1", 
    "AC106786.1", "RNF168", "CASP9", "MLXIPL", "CDX2", "ANKRD54", 
    "AC027682.7", "PRRG1", "AC209007.1", "AC026356.1", "AC104964.3", 
    "THG1L", "AC139530.1", "TEX22", "RPL7P59", "CRB2", "PDCD2", "HTR1D", 
    "PPIA", "YARS2", "AC092183.1", "TVP23C", "EEF1A1P6", "DNAAF1", 
    "INAFM1", "TRIOBP", "PIR", "NOP10", "DDX50P2", "ACTL8", "UQCRC2", 
    "RAB7A", "FBXL19-AS1", "PA2G4P4", "MAP4K2", "HADHAP2", "CKAP5", 
    "ZNF215", "GPX4", "MTX1", "TTLL4", "AHNAK", "KRT18P31", "ELP4", 
    "NCAM1", "DOK4", "ACAA2", "AC092431.2", "RPL22P8", "SP9", "ST6GALNAC6", 
    "TRA2A", "NOP2", "PLXNB2", "GADD45GIP1", "ASH1L-AS1", "ZNF843", 
    "RPS23P8", "KARS", "ZNF407", "TRAPPC2B", "ZNF518A", "KCNA3", 
    "UCKL1", "NPM1P6", "NOL3", "SPR", "NOP14", "LGMN", "VWA5A", "GAPDHP25", 
    "ARTN", "KBTBD7", "DCXR", "ANO6", "AC008147.4", "HAX1", "NEDD1", 
    "ZNF772", "AC020910.5", "DEPDC7", "HERC2P2", "DCHS1", "RPL37P1", 
    "PLIN1", "APLF", "ASB13", "NOP56", "PCSK4", "AMACR", "PLCL1", 
    "G0S2", "VASN", "PPP3CB", "PSMD4", "AL512791.2", "STK17A", "TSSK3", 
    "AC006460.2", "PPP1R26", "ELF1", "ACTBP8", "CSNK2B", "EDRF1", 
    "AC027097.1", "PTBP1", "KLLN", "GSTK1", "VKORC1L1", "COL11A2", 
    "C15orf62", "PRKCG", "ACVR2B-AS1", "RDH16", "WASHC2C", "AC114980.1", 
    "TXN2", "RPL19P16", "SRR", "RPS19BP1", "AC010680.5", "RN7SL832P", 
    "PCNA", "AC013403.2", "UBE2Q2", "AC016737.1", "TMPO", "API5P2", 
    "FAAH2", "STAMBP", "AL445487.1", "RPS2P44", "SRRM2-AS1", "AC079193.2", 
    "GATAD2A", "FEZF1", "HHEX", "RPF1", "IRX3", "AC093788.1", "STON2", 
    "TCN2", "CMSS1", "RPL27AP5", "DPH1", "AL136116.3", "TOMM5", "RPL14", 
    "AL365295.1", "AC090015.1", "VPS13B", "MXD4", "CHP2", "RN7SL535P", 
    "ACTG1P14", "AC087343.1", "PRPF38AP2", "GRIPAP1", "AC107068.1", 
    "CTSL", "CTSF", "SNHG11", "HHIP-AS1", "NMUR1", "R3HCC1L", "PPP1R13L", 
    "UGDH-AS1", "PDF", "GINS4", "PIPSL", "AC109347.1", "SLC16A10", 
    "RPS18P13", "ASPHD2", "EEF1E1P1", "AL096803.2", "RIC8B", "TFDP1", 
    "RBM6", "METTL4", "SULT4A1", "RPS2P46", "DSTN", "AC005072.1", 
    "ACAA1", "AC012313.1", "GSDME", "HNRNPA1P54", "TOGARAM2", "AC022498.2", 
    "TVP23B", "RPS12P26", "AC131235.1", "RCCD1", "INPP1", "WDFY1", 
    "AC026367.2", "MCF2L-AS1", "SETP21", "TMEM138", "SNORA33", "PABPN1", 
    "OGDH", "AL122023.1", "CHEK1", "CCN2", "AL139260.1", "SYMPK", 
    "CDKN2AIPNL", "MBD5", "NCOA4P4", "BLVRB", "DSE", "MIEF2", "RBM8B", 
    "ZNF350", "ADSL", "PARN", "LNCTAM34A", "PNMA6A", "PIGW", "TIMM29", 
    "IL21R-AS1", "NBPF3", "CYB5RL", "CERT1", "MAP7D2", "RUNX1", "AC127024.8", 
    "JAKMIP2", "NIBAN2", "TIGAR", "RNF227", "UNG", "MTOR", "TMEM107", 
    "JUND", "GEMIN7", "CALCOCO2", "TOR2A", "ATAD3A", "C1orf43", "CRK", 
    "AL353151.1", "GAPDHP23", "DEGS1", "KXD1", "AF131216.4", "GTF2H2C", 
    "LAMP2", "IPO7P1", "RAP2C-AS1", "LFNG", "MTND2P9", "CNOT7", "AC134349.1", 
    "GMPSP1", "NAA25", "MMP24OS", "VMP1", "NTNG2", "MCCC1", "RPS25", 
    "ILF2P2", "AL035530.2", "AC024451.1", "GNMT", "XRCC6P1", "NSMCE3", 
    "AC026124.1", "AC074194.1", "NHP2P1", "CHD2", "AC073046.1")
}

.exampleGeneset <- function(){
  head(.exampleBackground(),60)
}

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



# Loading RBP Position files in the correct format
.loadRBPPos <- function(species = c("human","mouse","rat")) {
  species <- match.arg(species)
  message("Only mouse at the moment, will be loaded at every species")
  RBPPos <- switch( species,
                    human = data("RBPPos_test"),
                    mouse = data("RBPPos_test"),
                    rat = data("RBPPos_test"), 
                    stop("No matched species"))
  
  colnames(RBPPos)[colnames(RBPPos)=="motif_id"] <- "Motif_ID"
  colnames(RBPPos)[colnames(RBPPos)=="sequence_name"] <- "Transcript ID"
  RBP_Names <- switch( species,
                       human = data("RBP_Names_Mus_Sample_CISBP"),
                       mouse = data("RBP_Names_Mus_Sample_CISBP"),
                       rat = data("RBP_Names_Mus_Sample_CISBP"), 
                       stop("No matched species"))
  RBPPos <- merge(RBPPos,RBP_Names_mus,by = "Motif_ID",all.x = TRUE)
  RBPPos
  
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

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

#' recapitalizeMiRs
#'
#' A utility function to reformat miRNA names so that there are no capital letters except the 'miR-'.
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
    x <- sign(dea$logFC)*-log10(dea$FDR)
  }else{
    x <- dea[[field]]
  }
  names(x) <- row.names(dea)
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
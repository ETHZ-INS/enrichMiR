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

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


.filterFamilies <- function(miRNAs, families){
  f <- families[intersect(miRNAs,names(families))]
  if(length(f)==0){
    names(families) <- sapply(names(families),FUN=function(x){ paste(strsplit(x,"-",fixed=T)[[1]][-1],collapse="-") })
    miRNAs <- sapply(miRNAs,FUN=function(x){ paste(strsplit(x,"-",fixed=T)[[1]][-1],collapse="-") })
    f <- families[intersect(miRNAs,names(families))]
  }
  return(f)
}


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

RNAfold <- function(seqs, bin.path="RNAfold", nthreads=1){
  f <- tempfile()
  write(as.character(seqs), f)
  a <- system(paste0('cat ', f, ' | ', bin.path, ' --jobs=',as.integer(nthreads),' --noPS'),
              intern=TRUE)
  a <- as.data.frame(t(vapply(strsplit(a[1:(length(a)/2)*2], " ( ",fixed=TRUE),
                FUN.VALUE=character(2L), FUN=identity)), stringsAsFactor=FALSE)
  colnames(a) <- c("fold","energy")
  a$energy <- as.numeric(gsub(")","",a$energy))
  a
}

#' optimize.ag
#'
#' Optimizes the a_g parameter based on matches' log_kd and a DEA
#'
#' @param matches 
#' @param dea 
#' @param ag.range 
#' @param costfn cost function; if null will use multiple ones
#' @param maximum Logical; whether to maximize
#'
#' @return if costfn is given, a list returning the optimized parameter and its objective;
#' otherwise the optimized parameter values based on the different cost functions.
#' @export
optimize.ag <- function(matches, dea, ag.range=10^c(-0.5,-4), costfn=NULL, maximum=TRUE){
  m2 <- matches[seqnames(matches) %in% row.names(dea),]
  if(is.null(costfn)){
    return(sapply(c("pearson","spearman","iMAD","NMI"), FUN=function(x){
      optimize.ag(matches=m2, dea=dea, costfn=x)$maximum
    }))
  }
  txs <- droplevels(factor(seqnames(m2)))
  lfc <- dea[levels(txs),"logFC.SVA.shrinked"]  # generalize this or change to something better
  if(!is.function(costfn))
    costfn <- switch(costfn,
                     "pearson"=function(x,y) cor(x,y,use="pairwise"),
                     "spearman"=function(x,y) cor(x,y,use="pairwise",method="spearman"),
                     "iMAD"=function(x,y) 1/median(abs(x-y),na.rm=TRUE),
                     "NMI"=function(x,y) aricode::NMI(x,y),
                     stop("Unknown cost function!"))
  f <- function(ag){
    x <- rowsum(ag/(ag-10^(m2$log_kd/-100)), txs)
    as.numeric(costfn(x[,1],lfc))
  }
  optimise(f, ag.range, maximum = maximum)
}

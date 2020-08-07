#' testEnrichment
#'
#' Creates an enrich.results object and performs enrichment analysis
#'
#' @param x The signature in which to look for a signal. This can be either:
#' \itemize{
#' \item A data.frame of the results of a differential expression analysis, 
#' with features (e.g. genes) as row names and with at least the following 
#' columns: `logFC`, `FDR`;
#' \item A numeric vector of the signal with feature (e.g. genes) as names, in
#' which case only tests based on a continuous signal will be possible;
#' \item A logical vector of membership to the geneset of interest, with feature
#' (e.g. genes) as names, in which case only tests based on a binary signal will
#' be available.
#' \item A vector of feature names; in this case the \code{background} argument will
#' also be required, and only tests based on a binary signal will be available.
#' }
#' @param sets Either a named list (or SimpleList) of features, or a data.frame
#' (or DataFrame) with at least the columns `set` and `feature` (if the columns 
#' `score` and `sites` are also present, they will be used by the 
#' appropriate tests).
#' @param background A character vector of background; ignored if `x` is not a
#' character vector.
#' @param tests Character vector of the tests to perform. See 
#' \link{\code{availableTests}} for the options.
#' @param sets.properties Any further information about the sets; this can 
#' either be a data.frame (or DataFrame), with row.names corresponding to names
#' of `sets` (or to alternative names), or a named vector (e.g. miRNA expression
#'  values). If this is given, sets not included within this object will be 
#'  discarded.
#' @param th.abs.logFC The minimum absolute log2-foldchange threshold for a 
#' feature to be considered differentially-expressed (default 0). Ignored if `x`
#' is not a DEA data.frame.
#' @param th.FDR The maximum FDR for a feature to be considered differentially-
#' expressed (default 0.05). Ignored if `x` is not a DEA data.frame.
#' @param minSize The minimum size of a set to be tested (default 5).
#' @param maxSize The maximum size of a set to be tested (default Inf).
#' @param gsea.maxSize The maximum number of targets for a miRNA family to be tested using GSEA (default 300).
#' @param gsea.permutations The number of permutations for GSEA (default 2000). See `?fgsea` for more information.
#' @param testOnlyAnnotated Whether to excluded features that are no set from 
#' the background (default FALSE).
#'
#' @return an enrich.results object.
#'
#' @export
testEnrichment <- function( x, sets, background=NULL, tests=NULL, 
                            sets.properties=NULL, th.abs.logFC=0, th.FDR=0.05,
                            minSize=5, maxSize=Inf, gsea.maxSize=1000, 
                            gsea.permutations=2000, testOnlyAnnotated=FALSE, 
                            keepAnnotation=FALSE, BPPARAM=NULL, 
                            field=NULL, ...){
  if(is.character(x)){
    if(is.null(background)) 
      stop("If `x` is a character vector, `background` should be given.")
    x <- background %in% x
    names(x) <- background
  }else{
    if(!is.null(background)) warning("`background` ignored.")
  }
  sets <- .list2DF(sets)
  
  if(is.null(dim(x))){
    if(testOnlyAnnotated) x <- x[names(x) %in% unique(sets$feature)]
    U <- names(x)
  }else{
    if(testOnlyAnnotated) x <- x[row.names(x) %in% unique(sets$feature),]
    U <- row.names(x)
    if(is.null(field)){
      if(!all(c("logFC","FDR") %in% colnames(x)))
        stop("`x` does not contain a 'logFC' and 'FDR' columns. Either rename ",
             "the columns or specify the `field` argument.")
    }else{
      if(!("field" %in% colnames(x)))
        stop("The `field` selected is not available in `x`.")
    }
  }
  if(is.null(U)) stop("`x` should be named.")

  atests <- availableTests(x, sets)
  if(is.null(tests)){
    tests <- setdiff(atests, c("gsea", "ks", "ks2", "mw"))
  }else{
    if(length(bad.tests <- setdiff(tolower(tests), atests))>0)
      stop("The following requested tests are not available with this kind of input: ", 
           paste(bad.tests, collapse=", "))
  }
  sets.properties <- .filterMatchSets(sets, sets.properties)
  # restrict sets according to properties and sizes
  sets <- sets[sets$set %in% row.names(sets.properties),]
  sets.properties$size <- table(sets$set)[row.names(sets.properties)]
  sets.properties <- sets.properties[sets.properties$size >= minSize & 
                                       sets.properties$size <= maxSize,
                                     , drop=FALSE]
  sets <- sets[sets$set %in% row.names(sets.properties),]
  if(is.factor(sets$set)) sets$set <- droplevels(sets$set)

  binary.signatures <- list()
  signature <- NULL
  if(!is.null(dim(x))){
    binary.signatures <- list(
      down=.dea2binary(x, th=th.FDR, th.alfc=th.abs.logFC, restrictSign=-1, ...),
      up=.dea2binary(x, th=th.FDR, th.alfc=th.abs.logFC, restrictSign=1, ...)
    )
    signature <- .dea2sig(x, field=field)
  }else if(is.logical(x)){
    binary.signatures <- list(x)
  }else if(is.numeric(x)){
    signature <- x
  }
  signature <- signature[!is.na(signature)]
  if(length(w <- which(is.infinite(signature)))>0)
    signature[w] <- max(abs(signature[-w]))*sign(signature[w])
  
  o <- new("enrich.results", input=list(x=x, sets.properties=sets.properties), 
           binary.signatures=binary.signatures, info=list(call=match.call()))
  o@overlaps <- lapply(binary.signatures, FUN=function(x){
    lapply(split(sets$feature, sets$set), y=names(x)[x], FUN=intersect)
  })
  if(keepAnnotation) o@input$sets <- sets

  message("Running the following tests: ", paste( tests, collapse=", "))
  names(tests2) <- tests2 <- setdiff(tests, "gsea")
  if(is.null(BPPARAM)){
    o@res <- unlist( lapply(tests2, sig=signature, sets=sets, 
                              binary.signatures=binary.signatures, 
                              FUN=.dispatchTest),
                     recursive=FALSE)
  }else{
    o@res <- unlist( bplapply(tests2, sig=signature, sets=sets, 
                            binary.signatures=binary.signatures, 
                            BPPARAM=BPPARAM, FUN=.dispatchTest),
                   recursive=FALSE)
  }
  if("gsea" %in% tests)
    o@res$gsea <- gsea(signature, sets, maxSize=gsea.maxSize, 
                       nperm=gsea.permutations, 
                       BPPARAM=ifelse(is.null(BPPARAM), SerialParam(), BPPARAM))
  
  return(o)
}


enrichMiR <- function( DEA, TS, miRNA.expression=NULL, families=NULL, cleanNames=FALSE, ...){
  if(is.null(families)){
    data("miR_families")
    families <- miR_families
  }
  if(cleanNames) names(families) <- sapply(names(families),FUN=.cleanMiRname)
  if(!is.null(miRNA.expression)){
    if(is.matrix(miRNA.expression) | is.data.frame(miRNA.expression)) miRNA.expression <- rowMeans(miRNA.expression,na.rm=T)
    if(cleanNames) names(miRNA.expression) <- sapply(names(miRNA.expression),FUN=.cleanMiRname)
    miRNA.expression <- miRNA.expression[which(miRNA.expression>0)]
  }
  colnames(TS) <- gsub("^family$","set", colnames(TS))
  TS <- DataFrame(TS)
  metadata(TS)$alternativeNames <- families
  testEnrichment(DEA, TS, sets.properties=miRNA.expression, ...)
}

#' getResults
#'
#' Returns the results table of an enrich.results object.
#'
#' @param object An object of class `enrichMiR`, as produced by the enrichMiR function.
#' @param test The name of a test, or a vector of such names. If ommitted (default), all available tests are returned. Note that more detailed results are returned when looking at a specific test.
#' @param nameCleanFun A function to clean the name of miRNAs. By default, nothing is done; you can try `nameCleanFun=.cleanMiRname` to remove prefixes.
#' 
#' @return a data.frame.
#'
#' @export
getResults <- function(object, test=NULL, nameCleanFun=function(x){x}){
  if(!is(object,"enrichMiR")) stop("object should be of class `enrichMiR'")
  if(!is.null(test) && !all(test %in% names(object@res))) 
    stop(paste("Required test not in the object's results. Available tests are:",paste(names(object@res),collapse=", ")))
  if(is.null(test) || length(test)==0) test <- names(object@res)
  ll <- object@res[test]
  if(length(ll)>1){
    ll <- lapply(ll,FUN=function(x){
      x[,grep("pvalue|FDR|enrichment|nes|beta|overlap",colnames(x)),drop=FALSE]
    })
    res <- .plMergeList(ll,all=T)
    rm(ll)
  }else{
    res <- ll[[1]]
  }

  ff <- split(names(object@families),object@families)
  ff <- sapply(ff,FUN=function(x){ paste(sapply(x,nameCleanFun),collapse=", ") })
  ff <- as.data.frame(ff,row.names=names(ff))
  colnames(ff) <- "miRNAs"
  ff <- ff[row.names(res),,drop=F]

  if(!is.null(object@miRNA.expression[[1]])){
    ff$expression <- object@miRNA.expression$family[row.names(ff)]
  }
  res <- cbind(ff,res)
  fdr <- grep("FDR",colnames(res))
  if(length(fdr)==1){
    fdr <- res[[fdr]]
  }else{
    fdr <- rowMeans(log10(as.matrix(res[,fdr])),na.rm=T)
  }
  res[order(fdr),]
}


availableTests <- function(x=NULL, sets=NULL){
  sigBinary <- sigContinuous <- setsScore <- TRUE
  if(!is.null(x)){
    sigBinary <- !is.null(dim(x)) || is.logical(x) || is.character(x)
    sigContinuous <- !is.null(dim(x)) || is.numeric(x)
  }
  if(!is.null(sets)){
    setsScore <- !is.null(dim(sets)) && "score" %in% colnames(sets)
    setsSites <- !is.null(dim(sets)) && "sites" %in% colnames(sets)
  }
  tests <- c("regmirb")
  if(sigBinary) tests <- c(tests, c("overlap","siteoverlap","woverlap"))
  if(sigContinuous) tests <- c(tests, c("mw","ks","ks2","gsea","areamir"))
  if(sigContinuous && setsScore) tests <- c(tests, c("modscore"))
  if(sigContinuous && setsSites) tests <- c(tests, c("modsites"))
  if(setsScore) tests <- c(tests, c("regmir"))
  sort(tests)
}

.filterMatchSets <- function(sets, props, tryPrefixRemoval){
  if(is.null(props)) return(data.frame(row.names=unique(sets$set)))
  if(is.null(dim(props))){
    if(is.null(names(props))){
      if(!is.character(props)) stop("Unrecognized set properties")
      props <- data.frame(row.names=props)
    }else{
      props <- as.data.frame(props)
    }
  }
  if(all(row.names(props) %in% sets$set)) return(props)
  
  if(isS4(sets) && !is.null(an <- metadata(sets)$alternativeNames) &&
     length(w <- which(row.names(props) %in% names(an)))>0){
    props <- aggregate(rbind(props[w,,drop=FALSE],props[-w,,drop=FALSE]),
                       by=list(set=c(an[row.names(props)[w]],row.names(props)[-w])),
                       FUN=function(x){
                         if(length(x)==1) return(x)
                         if(is.numeric(x)) return(x[which.max(abs(x),na.rm=TRUE)])
                         paste(sort(x), collapse=";")
                       })
    row.names(props) <- props[,1]
    props[,1] <- NULL
  }
  
  if(isS4(sets) && !any(row.names(props) %in% sets$set) && tryPrefixRemoval){
    names(an) <- sapply(strsplit(names(an),"-",fixed=T),FUN=function(x){ paste(x[-1],collapse="-") })
    metadata(sets)$alternativeNames <- an
    row.names(props) <- sapply(strsplit(row.names(props),"-",fixed=T),FUN=function(x){ paste(x[-1],collapse="-") })
    return(.filterMatchSets(sets, props, tryPrefixRemoval=FALSE))
  }

  if(!any(row.names(props) %in% sets$set))
    stop("The set properties do not correspond to the sets!")

  warning(sum(!(row.names(props) %in% sets$set)), " sets in the set.properties object",
          " are not in the sets...")
  props[intersect(row.names(props),sets$set),]
}

.list2DF <- function(sets){
  if(is(sets, "DataFrame")) return(sets)
  if(is.data.frame(sets)) return(DataFrame(sets))
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

.dispatchTest <- function(test, sig, sets, binary.signatures=list(), ...){
  ll <- list()
  if(length(binary.signatures)>0 && 
     test %in% c("overlap","siteoverlap","woverlap","regmir","regmirb")){
    ll <- lapply(binary.signatures, sets=sets, FUN=test)
  }
  if(!is.null(sig) && !(test %in% c("overlap","siteoverlap","woverlap"))){
    if(!is.null(dim(sig))) sig <- .dea2sig(sig, ...)
    ll[["continuous"]] <- get(test)(sig, sets)
  }
  ll
}
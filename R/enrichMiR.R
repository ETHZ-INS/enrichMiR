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
#' @param gsea.maxSize The maximum number of targets for a miRNA family to be 
#' tested using GSEA (default 300).
#' @param gsea.permutations The number of permutations for GSEA (default 2000).
#' See `?fgsea` for more information.
#' @param testOnlyAnnotated Whether to excluded features that are no set from 
#' the background (default FALSE).
#' @param field The field of `x` to use as continuous signal. If omitted, 
#' defaults to `sign(x$logFC)*-log10(x$FDR)`. Ignored if `x` is
#' not a data.frame.
#' @param BPPARAM \link{BiocParallel} multithreading parameters. Used to
#' multithread tests, and also within the `gsea` test.
#'
#' @return an enrich.results object.
#'
#' @import S4Vectors
#' @importFrom BiocParallel bplapply SerialParam
#' @export
testEnrichment <- function( x, sets, background=NULL, tests=NULL, 
                            sets.properties=NULL, th.abs.logFC=0, th.FDR=0.05,
                            minSize=5, maxSize=Inf, gsea.maxSize=1000, 
                            gsea.permutations=2000, testOnlyAnnotated=FALSE, 
                            keepAnnotation=FALSE, field=NULL, BPPARAM=NULL, 
                            ...){
  if(is.character(x)){
    if(is.null(background)) 
      stop("If `x` is a character vector, `background` should be given.")
    x <- background %in% x
    names(x) <- background
  }else{
    if(!is.null(background)) warning("`background` ignored.")
  }
  if(!is.null(dim(x))) x <- .homogenizeDEA(x)
  x <- .applySynonyms(x, sets)

  if( (is.null(dim(x)) && !any(names(x) %in% sets$feature)) ||
      (!is.null(dim(x)) && !any(row.names(x) %in% sets$feature)) )
    stop("There is no match between the features of the `sets` and those of `x`. ",
         "Are you sure that you are using the right annotation for your data?")
  
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
      stop("The following tests are not available with this kind of input: ", 
           paste(bad.tests, collapse=", "))
  }
  sets.properties <- .filterMatchSets(sets, sets.properties)
  # restrict sets according to properties and sizes
  sets <- sets[sets$set %in% row.names(sets.properties),]
  sets.properties$set_size <- table(sets$set)[row.names(sets.properties)]
  sets.properties <- sets.properties[sets.properties$set_size >= minSize & 
                                       sets.properties$set_size <= maxSize,
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
    FactorList(lapply(split(sets$feature, sets$set), 
                      FUN=function(y) sort(intersect(y,names(x)[x]))))
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


#' enrichMiR
#' 
#' A miRNA wrapper around `testEnrichment`, for continuity with previous 
#' enrichMiR versions.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features
#'  (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: 
#' `family`, `rep.miRNA`, `feature`, `sites`. Alternatiely, the output of an 
#' aggregated scan of miRNA kdModels, with columns `transcript`, `seed`, 
#' `score` (or `log_kd`), `n.8mer`, etc.
#' @param miRNA.expression A named vector of miRNAs expression values. 
#' miRNAs not in this vector are assumed to be not expressed in the system, and are not tested.
#' @param families A named vector of miRNA families, with individual miRNAs as names. 
#' If not given, internal data from the package will be used (mouse miRNA families from targetScan).
#' @param cleanNames Logical; whether to remove prefix from all miRNA names (default FALSE).
#' @param ... Passed to `testEnrichment`
#'
#' @return an enrich.results object.
#'
#' @export
#'
#' @examples
enrichMiR <- function( DEA, TS, miRNA.expression=NULL, families=NULL, cleanNames=FALSE, ...){
  if(is.null(families)){
    data("miR_families")
    families <- miR_families
  }
  if(cleanNames) names(families) <- sapply(names(families),FUN=.cleanMiRname)
  if(!is.null(miRNA.expression)){
    if(is.matrix(miRNA.expression) | is.data.frame(miRNA.expression))
      miRNA.expression <- rowMeans(miRNA.expression,na.rm=T)
    if(cleanNames) names(miRNA.expression) <- sapply(names(miRNA.expression),FUN=.cleanMiRname)
    miRNA.expression <- miRNA.expression[which(miRNA.expression>0)]
  }
  if(all(c("family", "feature") %in% colnames(TS))){
    colnames(TS) <- gsub("^family$","set", colnames(TS))
    TS <- DataFrame(TS)
    metadata(TS)$alternativeNames <- families
  }else if(all(c("seed","transcript","n.8mer") %in% colnames(TS))){
    TS$sites <- rowSums(TS[,grep("^n\\.", colnames(TS)),drop=FALSE], na.rm=TRUE)
    scoreField <- intersect(c("score", "repression", "log_kd"), colnames(TS))[1]
    TS <- TS[,c("seed","transcript","sites", scoreField)]
    colnames(TS) <- c("set","feature","sites","score")
    if(all(head(TS$score)>0)) TS$score <- -TS$score
  }
  TS <- DataFrame(TS)
  metadata(TS)$alternativeNames <- families
  testEnrichment(DEA, TS, sets.properties=miRNA.expression, ...)
}

#' getResults
#'
#' Returns the results table of an enrich.results object.
#'
#' @param object An object of class `enrich.results`.
#' @param test The name of a test, or a vector of such names. If ommitted (default), a 
#' merge of all available tests is returned. Note that more detailed results are returned
#'  when looking at a specific test.
#' @param getFeatures Logical; whether to return overlapping features.
#' @param flatten Logical; whether to return a normal data.frame instead of a DataFrame.
#' 
#' @return a DataFrame or data.frame.
#'
#' @export
getResults <- function(object, test=NULL, getFeatures=TRUE, flatten=FALSE){
  if(!is(object,"enrich.results")) stop("object should be of class `enrich.results'")
  if(is.null(test) && length(object@res)==1) test <- names(object@res)
  if(!is.null(test) && !all(test %in% names(object@res))) 
    stop(paste("Required test not in the object's results. Available tests are:",
               paste(names(object@res),collapse=", ")))
  if(is.null(test) || length(test)==0) test <- names(object@res)
  ll <- object@res[test]
  if(length(ll)>1){
    ll <- lapply(ll,FUN=function(x){
      x[,grep("pvalue|FDR|enrichment|nes|beta|coefficient|overlap",colnames(x)),drop=FALSE]
    })
    res <- .plMergeList(ll,all=T)
    rm(ll)
  }else{
    res <- ll[[1]]
  }

  res <- cbind(object@input$sets.properties[row.names(res),,drop=FALSE], res)

  fdr <- grep("FDR",colnames(res))
  if(length(fdr)==1){
    ob <- res[[fdr]]
  }else{
    ob <- 10^rowMeans(log10(10^-16+as.matrix(res[,fdr])),na.rm=T)
    res$FDR.mean <- rowMeans(as.matrix(res[,fdr]),na.rm=T)
    res$FDR.geomean <- ob
  }
  res <- dround(res[order(ob),], roundGreaterThan1=TRUE)
  if(getFeatures && length(object@overlaps)>0){
    for(f in names(object@overlaps)){
      if(flatten){
        res[[paste0("genes.", f)]] <- sapply(object@overlaps[[f]][row.names(res)], collapse=",", FUN=paste)
      }else{
        res <- DataFrame(res)
        res[[paste0("genes.", f)]] <- FactorList(lapply(row.names(res), FUN=function(x) object@overlaps[[f]][[x]]))
      }
    }
  }
  if(flatten) {
    return(res)
  }else{
    res <- DataFrame(res) 
    }
  res
}


#' availableTests
#'
#' Returns the available enrichment tests, for the given inputs (if specified)
#'
#' @param x Optional `x` input for `testEnrichment`
#' @param sets Optional `sets` input for `testEnrichment`
#'
#' @return A character vector of available tests
#' @export
#'
#' @examples
#' availableTests()
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

.filterMatchSets <- function(sets, props, tryPrefixRemoval=TRUE){
  if(is.null(props)){
    if(is.null(fams <- metadata(sets)$families))
      return(data.frame(row.names=unique(sets$set)))
    props <- split(names(fams),fams)
    props <- data.frame(row.names=names(props),
                        members=sapply(props, FUN=function(x) paste(sort(x),collapse=";")))
    if(length(miss <- setdiff(unique(sets$set),row.names(props)))>0){
      props <- rbind(props, data.frame(row.names = miss, members=miss))
    }
    return(props)
  }
  if(is.null(dim(props))){
    if(is.null(names(props))){
      if(!is.character(props)) stop("Unrecognized set properties")
      props <- data.frame(row.names=props)
    }else{
      props <- as.data.frame(props)
    }
  }
  usets <- unique(sets$set)
  if(all(row.names(props) %in% usets)) return(props)
  if(isS4(sets)){
    if(is.null(metadata(sets)$alternativeNames) && !is.null(metadata(sets)$families))
      metadata(sets)$alternativeNames <- metadata(sets)$families
    if(!is.null(an <- metadata(sets)$alternativeNames)){
      props <- .agDF(props, an, TRUE)
      if(all(row.names(props) %in% usets)) return(props)
    }
    if(!is.null(an) && tryPrefixRemoval &&
       !any(row.names(props) %in% usets) ){
      names(an) <- sapply(strsplit(names(an),"-",fixed=TRUE),FUN=function(x){ paste(x[-1],collapse="-") })
      metadata(sets)$alternativeNames <- an
      props <- .agDF(props, sapply(strsplit(row.names(props),"-",fixed=T),FUN=function(x) paste(x[-1],collapse="-") ), FALSE)
      return(.filterMatchSets(sets, props, tryPrefixRemoval=FALSE))
    }
  }
  if(!any(row.names(props) %in% usets))
    stop("The set properties do not correspond to the sets!")
  warning(sum(!(row.names(props) %in% usets)), " sets in the set.properties object",
          " are not in the sets...")
  props[intersect(row.names(props),usets),,drop=FALSE]
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
    res <- get(test)(sig, sets)
    if(test %in% c("regmir","regmirb")){
      ll[["continuous"]] <- res
    }else{
      ll <- c(ll, list(res))
    }
    ll
  }
  ll
}

.applySynonyms <- function(x, sets){
  if(is.null(metadata(sets)) || is.null(metadata(sets)$feature.synonyms)) return(x)
  ss <- metadata(sets)$feature.synonyms
  if(is.null(dim(x)) && any(names(x) %in% names(ss))){
    n <- names(x)
    w <- which(n %in% names(ss))
    n[w] <- as.character(ss[n[w]])
    if((nDup <- sum(duplicated(n)))>0){
      warning(nDup, " duplicated feature(s) dropped.")
      w <- which(duplicated(n))
      x <- x[-w]
      names(x) <- n[-w]
    }else{
      names(x) <- n
    }
  }else if(!is.null(dim(x)) && any(row.names(x) %in% names(ss))){
    # x <- .homogenizeDEA(x)
    x <- x[order(x$FDR),]
    n <- row.names(x)
    w <- which(n %in% names(ss))
    n[w] <- as.character(ss[n[w]])
    w <- duplicated(n)
    x <- x[!w,]
    row.names(x) <- n[!w]
  }
  return(x)
}

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
#' \code{\link{enrichMiR::availableTests}} for the options.
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
#' @param doCheck Whether to check the inputs. This argument is for internal use
#' and should not be disabled!
#'
#' @return an enrich.results object.
#'
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom IRanges FactorList CharacterList
#' @importFrom S4Vectors metadata metadata<- DataFrame
#' @export
testEnrichment <- function( x, sets, background=NULL, tests=NULL, 
                            sets.properties=NULL, th.abs.logFC=0, th.FDR=0.05,
                            minSize=5, maxSize=Inf, gsea.maxSize=3000, 
                            gsea.permutations=4000, testOnlyAnnotated=FALSE, 
                            keepAnnotation=FALSE, field=NULL, BPPARAM=NULL, 
                            familyRound=NULL, checkSynonyms=TRUE, doCheck=TRUE, 
                            ...){
  if(doCheck){
    if(is.character(x)){
      if(is.null(background)) 
        stop("If `x` is a character vector, `background` should be given.")
      x <- background %in% x
      names(x) <- background
    }else{
      if(!is.null(background)) warning("`background` ignored.")
    }
    message("Preparing inputs...")
    if(!is.null(dim(x))) x <- .homogenizeDEA(x)
    if(checkSynonyms) x <- .applySynonyms(x, sets)

    sets <- .filterInput(sets, x, min.size=minSize, max.size=maxSize)
    x <- sets$signal
    sets <- sets$sets
  }
  
  atests <- availableTests(x, sets)
  if(is.null(tests)){
    tests <- intersect(c("siteoverlap","areamir"), atests)
  }else{
    tests <- tolower(tests)
    if("regmir" %in% tests)
      tests <- c(setdiff(tests, "regmir"), grep("regmir",atests,value=TRUE))
    if(length(bad.tests <- setdiff(tests, atests))>0)
      stop("The following tests are not available with this kind of input: ", 
           paste(bad.tests, collapse=", "))
  }
  
  sets.properties <- .filterMatchSets(sets, sets.properties)
  # restrict sets according to properties and sizes
  sets <- sets[sets$set %in% row.names(sets.properties),]

  setsizes <- table(sets$set)
  sets.properties$set_size <- as.integer(setsizes[row.names(sets.properties)])
  if(!is.null(sets.properties$members) && 
     length(w <- is.na(sets.properties$set_size))>0)
    sets.properties$set_size[w] <- 
    unlist(lapply(strsplit(sets.properties$members[w], ";"), FUN=function(x){
      sum(as.integer(setsizes[x]),na.rm=TRUE)
    }))
  setsizes <- setsizes[setsizes>=minSize & setsizes<=maxSize]
  sets <- sets[sets$set %in% names(setsizes),]
  
  if("woverlap" %in% tests && length(setsizes)<10){
    tests <- setdiff(tests, "woverlap")
    message("Skipping 'woverlap' -- too few sets in target.")
  }
  
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
           binary.signatures=binary.signatures, info=list(call=match.call(),
                                                          type = "enrichment"))
  if(any(grep("-miR-|-mir-|-let-",head(row.names(sets.properties)))) && 
     any(tests %in% c("siteoverlap","overlap","siteoverlap2","woverlap")) && 
     !is.null(fams <- metadata(sets)$families)){
    o@input <- c(o@input,list(families = fams))
  }
  if(is.null(names(binary.signatures))) names(binary.signatures) <- "features"
    o@overlaps <- lapply(binary.signatures, FUN=function(x){
      FactorList(lapply(split(sets$feature, sets$set), 
                        FUN=function(y) sort(intersect(y,names(x)[x]))))
  })
  if(keepAnnotation) o@input$sets <- sets

  message("Running the following tests: ", paste( tests, collapse=", "))
  names(tests2) <- tests2 <- setdiff(tests, "gsea")
  if(length(tests2)>0){
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
  }
  if("gsea" %in% tests){
    if(is.null(BPPARAM)) BPPARAM <- SerialParam()
    o@res$gsea <- gsea(signature, sets, maxSize=gsea.maxSize, 
                       nperm=gsea.permutations,  BPPARAM=BPPARAM)
  }
  
  return(o)
}


#' enrichMiR
#' 
#' A miRNA wrapper around `testEnrichment`, for continuity with previous 
#' enrichMiR versions. It is recommended to use \code{\link{testEnrichment}} 
#' instead.
#'
#' @param DEA A data.frame of the results of a differential expression analysis,
#' with features (e.g. genes) as row names and with at least the following 
#' columns: `logFC`, `FDR`
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
enrichMiR <- function( DEA, TS, miRNA.expression=NULL, families=NULL, 
                       cleanNames=FALSE, ...){
  if(is.null(families)) families <- metadata(TS)$families
  if(cleanNames) names(families) <- sapply(names(families),FUN=.cleanMiRname)
  if(!is.null(miRNA.expression)){
    if(is.matrix(miRNA.expression) | is.data.frame(miRNA.expression))
      miRNA.expression <- rowMeans(miRNA.expression,na.rm=T)
    if(cleanNames) names(miRNA.expression) <- sapply(names(miRNA.expression),
                                                     FUN=.cleanMiRname)
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
    ff <- lapply(ll,FUN=function(x){
      x[,grep("pvalue|FDR|enrichment|nes|beta|coefficient|overlap",colnames(x)),drop=FALSE]
    })
    res <- .plMergeList(ff,all=T)
    rm(ff)
  }else{
    res <- ll[[1]]
  }
  if(object@info$type == "colocalization"){
    if(!all(c("miRNA","Partner") %in% colnames(res))){
      res <- cbind(ll[[1]][,c(1,2)],res)}
    res$miRNA_name <- object@input$sets.properties$miRNA[res$miRNA,"members"]
    if(length(object@input$sets.properties) == 1){
      res$Partner_name <- object@input$sets.properties$miRNA[res$Partner,"members"]
    }else{
      res$Partner_name <- object@input$sets.properties$Partner[res$Partner,"members"]
    }
  }else{
    if(any(grepl("-miR-|-mir-|-let-",head(row.names(object@input$sets.properties)))) && 
       any(test %in% c("siteoverlap.down","overlap.down","siteoverlap2.down",
                       "siteoverlap.up","siteoverlap2.up","siteoverlap.features",
                       "overlap.up","overlap.features","siteoverlap2.features")) && 
       !is.null(object@input$families)){
      fams <- object@input$families
      sets.properties <- object@input$sets.properties
      props1 <- data.frame(row.names = names(fams),seed = as.character(fams))
      sets.properties <- merge(sets.properties,props1, by = 0, all.x = TRUE)
      cols <- setdiff(colnames(sets.properties), c("Row.names", "set_size", "seed"))
      sets.properties <- setNames(aggregate(sets.properties[,cols], 
                                            by=list(sets.properties$seed), FUN=sum),
                                  c("seed",cols))
      props2 <- split(names(fams),fams)
      props2 <- data.frame(row.names=names(props2),
                           members=sapply(props2, FUN=function(x) paste(sort(x),collapse=";")))
      
      sets.properties <- merge(sets.properties,props2, by.x = "seed", by.y = 0, all.x = TRUE)
      row.names(sets.properties) <- sets.properties$seed
      sets.properties$seed <- NULL
      object@input$sets.properties <- sets.properties
      }
    res <- cbind(object@input$sets.properties[row.names(res),,drop=FALSE], res)
  }
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
        if(object@info$type == "colocalization"){
          res[[paste0(f,"_miRNA")]] <- sapply(object@overlaps[[f]][res$miRNA], collapse=",", FUN=paste)
        }else{
          if(all(row.names(res) %in% names(object@overlaps[[f]]))){
            res[[paste0("genes.", f)]] <- sapply(object@overlaps[[f]][row.names(res)], collapse=",", FUN=paste)
          }else if(!is.null(res$members)){
            # this is quite slow...
            ov <- sapply(strsplit(res$members, ";"), FUN=function(x){
              x <- object@overlaps[[f]][intersect(x, names(object@overlaps[[f]]))]
              paste(sort(unique(unlist(x, use.names=FALSE))), collapse=",")
            })
          }
        }
      }else{
        res <- DataFrame(res)
        if(object@info$type == "colocalization"){
          res[[paste0(f,"_miRNA")]] <- FactorList(lapply(res$miRNA, FUN=function(x) object@overlaps[[f]][[x]]))
        }else{
          res[[paste0("genes.", f)]] <- FactorList(lapply(row.names(res), FUN=function(x) object@overlaps[[f]][[x]]))
        }
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
  sigBinary <- sigContinuous <- setsScore <- setsSites <- TRUE
  atests <- NULL
  if(!is.null(x)){
    sigBinary <- !is.null(dim(x)) || is.logical(x) || is.character(x)
    sigContinuous <- !is.null(dim(x)) || is.numeric(x)
  }
  if(!is.null(sets)){
    if(.is.matrix(sets)){
      content <- attr(sets, "content")
      if(is.null(content) || !(content %in% c("sites","score"))){
        setsSites <- setsScore <- !is.logical(sets[1,1])  
      }else{
        setsSites <- content=="sites"
        setsScore <- content=="score"
      }
    }else{
      setsScore <- !is.null(dim(sets)) && "score" %in% colnames(sets)
      setsSites <- !is.null(dim(sets)) && "sites" %in% colnames(sets)
      if(is(sets, "DFrame")) atest <- metadata(sets)$tests
    }
  }
  if(!setsScore && !setsSites) return(c("overlap"))
  tests <- c("regmir.bb")
  if(sigBinary) tests <- c(tests, c("overlap","siteoverlap","woverlap"))
  if(sigContinuous) tests <- c(tests, c("mw","ks","gsea","areamir"))
  if(sigContinuous && setsScore) tests <- c(tests, c("modscore","regmir.cc","ebayes","lmadd"))
  if(sigContinuous && setsSites) tests <- c(tests, c("modsites"))
  if(setsScore && sigBinary) tests <- c(tests, c("regmir.bc"))
  if(!is.null(atests)) tests <- intersect(tests,atests)
  sort(tests)
}

.filterInput <- function(sets, signal, min.size=5, max.size=Inf, testOnlyAnnotated=FALSE){
  if(!is.null(dim(signal))){
    features <- row.names(signal)
  }else{
    features <- names(signal)
  }
  stopifnot(!is.null(features))
  if(.is.matrix(sets)){
    stopifnot(!is.null(row.names(sets)) && !is.null(colnames(sets)))
    if(testOnlyAnnotated)  features <- intersect(features, row.names(sets))
    sets <- sets[row.names(sets) %in% features,,drop=FALSE]
    cs <- colSums(sets!=0, na.rm=TRUE)
    sets[,cs>=min.size & cs<max.size,drop=FALSE]
    if(ncol(sets)==0) stop(.noMatchingFeatures())
    return(list(sets=sets, signal=signal[row.names(sets)]))
  }
  sets <- .list2DF(sets)
  sets$feature <- as.factor(sets$feature)
  sets$set <- as.factor(sets$set)
  if(testOnlyAnnotated) features <- intersect(levels(sets$feature), features)
  i <- which(levels(sets$feature) %in% features)
  sets <- sets[as.integer(sets$feature) %in% i,]
  sets$feature <- droplevels(sets$feature)
  tt <- table(sets$set)
  i <- which(levels(sets$set) %in% names(tt)[tt>=min.size & tt<max.size])
  if(length(i)==0) stop(.noMatchingFeatures())
  sets <- sets[as.integer(sets$set) %in% i,]
  sets$set <- droplevels(sets$set)
  if(!is.null(dim(signal))){
    return(list(sets=sets, signal=signal[row.names(signal) %in% levels(sets$feature),])) 
  }else{
    return(list(sets=sets, signal=signal[levels(sets$feature)])) 
  }
}

.noMatchingFeatures <- function(){
  stop("There is no match between the features of the `sets` and those of `x`. ",
       "Are you sure that you are using the right annotation for your data?")
}

.filterMatchSets <- function(sets, props, tryPrefixRemoval=TRUE){
  if(isS4(sets)){
    fams <- metadata(sets)$families
  }else{
    fams <- attr(sets, "families")
  }
  if(.is.matrix(sets))
    sets <- data.frame(feature=row.names(sets), set=colnames(sets))
  if(is.null(props)){
    if(is.null(fams)) return(data.frame(row.names=unique(sets$set)))
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
      props <- as.data.frame(props, row.names=names(props))
    }
  }
  usets <- unique(sets$set)
  if(all(row.names(props) %in% usets)) return(props)
  if(!any(grepl("-miR-|-mir-|-let-",head(usets)))){
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
  }
  if(!any(row.names(props) %in% usets))
    stop("The set properties do not correspond to the sets!")
  warning(sum(!(row.names(props) %in% usets)), " sets in the set.properties object",
          " are not in the sets...")
  props[intersect(row.names(props),usets),,drop=FALSE]
}


.dispatchTest <- function(test, sig, sets, binary.signatures=list(), ...){
  ll <- list()
  if(length(binary.signatures)>0 && 
     test %in% c("overlap","siteoverlap","woverlap","regmir.bb","regmir.bc")){
    ll <- lapply(binary.signatures, FUN=function(sig){
      tryCatch(get(test)(sig, sets=sets), 
               error=function(e){
                 message("Test '",test,"' failed with:")
                 message(e)
                 return(NULL)
               })
      })
  }
  if(!is.null(sig) && !(test %in% c("overlap","siteoverlap","woverlap"))){
    if(!is.null(dim(sig))) sig <- .dea2sig(sig, ...)
    res <- tryCatch(get(test)(sig, sets=sets), 
                     error=function(e){
                       message("Test '",test,"' failed with:")
                       message(e)
                       return(NULL)
                     })
    ll <- c(ll, list(res))
  }
  ll
}

.applySynonyms <- function(x, sets){
  if(is.null(metadata(sets)) || is.null(metadata(sets)$feature.synonyms))
    return(x)
  ss <- metadata(sets)$feature.synonyms
  if(is.null(dim(x)) && any(w <- grepl("^ENS",(names(x)))))
    names(x)[w] <- gsub("\\.[0-9]+$", "", names(x)[w])
  if(!is.null(dim(x)) && any(w <- grepl("^ENS",(row.names(x)))))
    row.names(x)[w] <- gsub("\\.[0-9]+$", "", row.names(x)[w])
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
    w <- duplicated(n) | is.na(n)
    x <- x[!w,]
    row.names(x) <- n[!w]
  }
  return(x)
}

#' enrichMiR
#'
#' Creates an enrichMiR object and performs miRNA target enrichment analyses
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features
#'  (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: 
#' `family`, `rep.miRNA`, `feature`, `sites`.
#' @param miRNA.expression A named vector of miRNAs expression values. 
#' miRNAs not in this vector are assumed to be not expressed in the system, and are not tested.
#' @param families A named vector of miRNA families, with individual miRNAs as names. 
#' If not given, internal data from the package will be used (mouse miRNA families from targetScan).
#' @param th.abs.logFC The minimum absolute log2-foldchange threshold for a feature/gene to be 
#' considered differentially-expressed (default 0).
#' @param th.FDR The maximum FDR for a feature/gene to be considered differentially-expressed (default 0.05).
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#' @param gsea.maxSize The maximum number of targets for a miRNA family to be tested using GSEA (default 300).
#' @param gsea.permutations The number of permutations for GSEA (default 2000). See `?fgsea` for more information.
#' @param gsea.fdr.thres The FDR threshold for inclusion in GSEA (default 0.2).
#' @param testOnlyAnnotated Whether to excluded features that are bound by no miRNA (default FALSE).
#' @param tests Character vector of the tests to perform. Any combination of: 
#' `overlap`, `wo` (weighted overlap), `siteMir`, `KS`, `KS2`, `MW`, `GSEA`, `modSites`, `modScore`, 
#' or NULL to perform all tests (default).
#' @param cleanNames Logical; whether to remove prefix from all miRNA names (default FALSE).
#'
#' @return an enrichMiR object.
#'
#' @export
enrichMiR <- function( DEA, TS, miRNA.expression=NULL, families=NULL, 
                       th.abs.logFC=0, th.FDR=0.05, minSize=5, gsea.maxSize=1000,
                       gsea.permutations=2000, gsea.fdr.thres=0.2, 
                       testOnlyAnnotated=FALSE, tests=NULL, cleanNames=FALSE, ...){
    if(!is.null(tests)) tests <- match.arg(tolower(tests), choices=c("overlap","sitemir","wo","mw","ks","ks2","gsea","modscore","modsites","areamir","areamirp","regmir","regmirb"), several.ok = T)
    if(is.null(families)){
        data("miR_families")
        families <- miR_families
        if(cleanNames) names(families) <- sapply(names(families),FUN=.cleanMiRname)
    }
    if(!is.null(miRNA.expression)){
        if(is.matrix(miRNA.expression) | is.data.frame(miRNA.expression)) miRNA.expression <- rowMeans(miRNA.expression,na.rm=T)
        if(cleanNames) names(miRNA.expression) <- sapply(names(miRNA.expression),FUN=.cleanMiRname)
        miRNA.expression <- miRNA.expression[which(miRNA.expression>0)]
        families <- .filterFamilies(names(miRNA.expression), families)
        tmp <- aggregate(miRNA.expression[names(families)],by=list(family=families),na.rm=T,FUN=sum)
        fam.expr <- tmp[,2]
        names(fam.expr) <- tmp[,1]
        miRNA.expression <- list(family=fam.expr, miRNA=miRNA.expression)
    }else{
        miRNA.expression <- list(family=NULL, miRNA=NULL)
    }
    down <- .dea2binary(DEA, th=th.FDR, th.alfc=th.abs.logFC, restrictSign=-1, ...)
    up <- .dea2binary(DEA, th=th.FDR, th.alfc=th.abs.logFC, restrictSign=1, ...)

    TS <- TS[which(as.character(TS$family) %in% families),]
    o <- new("enrichMiR", DEA=DEA, TS=as.data.frame(TS), families=families, miRNA.expression=miRNA.expression, info=list(call=match.call()))
    if(is.null(tests) || "areamir" %in% tests) o@res$aREAmir=aREAmir(DEA, TS, minSize, pleiotropy=FALSE)
    if(is.null(tests) || "areamirp" %in% tests) o@res$aREAmir2=aREAmir(DEA, TS, minSize, pleiotropy=TRUE)
    TS <- as.data.frame(TS)
    if(is.null(tests) || "overlap" %in% tests) o@res$EN.up <- EA(up, TS, minSize, testOnlyAnnotated)
    if(is.null(tests) || "overlap" %in% tests) o@res$EN.down <- EA(down, TS, minSize, testOnlyAnnotated)
    #if(is.null(tests) || "overlap" %in% tests) o@res$EN.combined <- .combTests(o@res$EN.up, o@res$EN.down)    
    if(is.null(tests) || "wo" %in% tests) o@res$wEN.up <- wEA(up, TS, minSize, testOnlyAnnotated)
    if(is.null(tests) || "wo" %in% tests) o@res$wEN.down <- wEA(down, TS, minSize, testOnlyAnnotated)
    #if(is.null(tests) || "wo" %in% tests) o@res$wEN.combined <- .combTests(o@res$wEN.up, o@res$wEN.down)
    if(is.null(tests) || "siteMir" %in% tests) o@res$siteMir.up <- siteMir(up, TS, minSize, testOnlyAnnotated)
    if(is.null(tests) || "siteMir" %in% tests) o@res$siteMir.down <- siteMir(down, TS, minSize, testOnlyAnnotated)
    #if(is.null(tests) || "siteMir" %in% tests) o@res$siteMir.combined <- .combTests(o@res$siteMir.up, o@res$siteMir.down)
    if(is.null(tests) || "mw" %in% tests) o@res$MW=MW(DEA, TS, minSize)
    if(is.null(tests) || "ks" %in% tests) o@res$KS=KS(DEA, TS, minSize)
    if(is.null(tests) || "ks2" %in% tests) o@res$KS2=KS2(DEA, TS, minSize)
    if(is.null(tests) || "gsea" %in% tests) o@res$GSEA=gsea(DEA, TS, minSize, maxSize=gsea.maxSize, nperm=gsea.permutations, fdr.thres=gsea.fdr.thres)
    if(is.null(tests) || "modsites" %in% tests) o@res$modSites=plMod(DEA, TS, minSize, var="sites", correctForLength=T)
    if(is.null(tests) || "modscore" %in% tests) o@res$modScore=plMod(DEA, TS, minSize, var="score", correctForLength=F)
    if(is.null(tests) || "regmir" %in% tests) o@res$regmir.down=regmir(down, TS, binary=FALSE, keepAll=TRUE)
    if(is.null(tests) || "regmir" %in% tests) o@res$regmir.up=regmir(up, TS, binary=FALSE, keepAll=TRUE)
    if(is.null(tests) || "regmirb" %in% tests) o@res$regmirb.down=regmir(down, TS, binary=TRUE, keepAll=TRUE)
    if(is.null(tests) || "regmirb" %in% tests) o@res$regmirb.up=regmir(up, TS, binary=TRUE, keepAll=TRUE)
    return(o)
}

.combTests <- function(up,down){
  library(aggregation)
  if( !identical(colnames(up),colnames(down)) || !identical(dim(up),dim(down)) ){
    stop("The tests should be of the same kind and have the same rows.")
  }
  ff <- intersect(c("overlap","enrichment","over.pvalue","under.pvalue","features"),colnames(up))
  m <- merge( down[,c("annotated",ff)], 
              up[,ff], by="row.names", all=T, suffixes=c(".down",".up"))
  row.names(m) <- m[,1]
  m$enrichment <- apply(m[,c("enrichment.down","enrichment.up")],1,FUN=function(x){ 
    x <- mean(c(-x[1],x[2]),na.rm=T)
    if(is.na(x)) return(0)
    x
  })
  m$overlap <- rowSums(as.matrix(m[,grep("overlap",colnames(m))]))
  m$comb.pvalue <- suppressWarnings(apply(m[,c("enrichment",colnames(m)[grep("pvalue",colnames(m))])],1,FUN=function(x){ 
    if(x["enrichment"]>0){
      x <- x[c("under.pvalue.down","over.pvalue.up")]
    }else{
      x <- x[c("over.pvalue.down","under.pvalue.up")]
    }
    fisher(x)
  }))
  m <- m[,grep("under.pvalue|over.pvalue",colnames(m),invert=T)]
  m$FDR <- p.adjust(m$comb.pvalue,method="fdr")
  m[order(m$FDR,m$comb.pvalue),-1]
}

#' enrichMiR.results
#'
#' Returns the results table of an enrichMiR object.
#'
#' @param object An object of class `enrichMiR`, as produced by the enrichMiR function.
#' @param test The name of a test, or a vector of such names. If ommitted (default), all available tests are returned. Note that more detailed results are returned when looking at a specific test.
#' @param nameCleanFun A function to clean the name of miRNAs. By default, nothing is done; you can try `nameCleanFun=.cleanMiRname` to remove prefixes.
#' 
#' @return a data.frame.
#'
#' @export
enrichMiR.results <- function(object, test=NULL, nameCleanFun=function(x){x}){
  if(!is(object,"enrichMiR")) stop("object should be of class `enrichMiR'")
  if(!is.null(test) && !all(test %in% names(object@res))) stop(paste("Required test not in the object's results. Available tests are:",paste(names(object$res),collapse=", ")))
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
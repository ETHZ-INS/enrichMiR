.plMergeList <- function(ll,...){
  if(!is.null(names(ll))) for(i in 1:length(ll)){
    colnames(ll[[i]]) <- paste(names(ll)[i],colnames(ll[[i]]),sep=".")
  }
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
#' @param x A character vector or an array (in which case the function will be 
#' applied to row.names)
#' @param gformat The gene format, either 'human' (all caps) or 'mouse' (first 
#' letter capitalized)
#' 
#' @return An object of the same type and dimensions as `x`.
#' 
#' @importFrom tools toTitleCase
#' @export
recapitalizeGenes <- function(x, gformat="mouse"){
  gformat <- match.arg(gformat, choices=c("human","mouse"))
  if(is(x,'data.frame') | is(x,'matrix')){
    row.names(x) <- recapitalizeGenes(row.names(x), gformat)
    return(x)
  }
  switch(gformat,
    mouse=tools::toTitleCase(tolower(x)),
    human=sapply(x,toupper)
  )
}

#' recapitalizeMiRs
#'
#' A utility function to reformat miRNA names so that there are no capital 
#' letters except the 'miR-'.
#'
#' @param x A character vector or an array (in which case the function will be 
#' applied to row.names)
#' 
#' @return An object of the same type and dimensions as `x`.
#'
#' @export
recapitalizeMiRs <- function(x){
  if(is(x,'data.frame') | is(x,'matrix')){
    if(is.null(row.names(x)))
      stop( "Given an array without row.names... please apply to the column ",
            "containing the miRNA names.")
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
  dea <- dea[order(dea$FDR),]
  if("PValue" %in% colnames(dea))
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
    if("PValue" %in% colnames(dea)){
      x <- sign(dea$logFC)*-log10(dea$PValue)
    }else{
      x <- sign(dea$logFC)*-log10(dea$FDR)
    }
  }else{
    dea <- dea[!is.na(dea[[field]]),]
    x <- dea[[field]]
  }
  names(x) <- row.names(dea)
  x
}

.homogenizeDEA <- function(x, keepTop=TRUE){
  if(is(x,"data.table")){
    if(any(duplicated(x[[1]]))){
      if(keepTop){
        x <- x[order(x[[head(grep("padj|adj\\.P\\.Val|q_value|qval|FDR", colnames(x)),1)]]),]
        x <- x[!duplicated(x[[1]]),]
      }else{
        x <- aggregate(x[,-1,drop=FALSE], by=list(gene=x[[1]]), FUN=mean)
      }
    }
    x <- x[!is.na(x[[1]]),]
    x <- data.frame(row.names=x[[1]], as.data.frame(x[,-1,drop=FALSE]))
  }
  x <- as.data.frame(x)
  if(length(w <- grep("^ENS",row.names(x)))>0)
    row.names(x)[w] <- gsub("\\..*","",row.names(x)[w])

  colnames(x) <- gsub("log2FoldChange|log2Fold|log2FC|log2\\(fold_change\\)|log2\\.fold_change\\.",
                      "logFC", colnames(x))
  
  abf <- colnames(df)[which(colnames(df) %in% c("meanExpr", "AveExpr", 
                                                "baseMean", "logCPM"))]
  if (length(abf) == 1) {
    x$meanExpr <- x[, abf]
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
  x <- x[!is.na(x$logFC) & !is.na(x$FDR),]
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
  if(is.null(names(new.rn))) return(df)
  tt <- split(names(new.rn),new.rn)
  tt <- data.frame(row.names=names(tt),
                   members=sapply(tt, FUN=function(x) paste(sort(x),collapse=";")))
  if(ncol(df)==0) return(data.frame(row.names=unique(new.rn)))
  if(!any(duplicated(new.rn))){
    row.names(df) <- new.rn
    df <- merge(df,tt,by = 0, all.x = TRUE)
    row.names(df) <- df$Row.names
    return(df[,-c(1),drop=FALSE])
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

#' @importFrom S4Vectors DataFrame metadata metadata<- 
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
  library(SummarizedExperiment)
  data("microRNAome", package="microRNAome")
  if(is.null(x)) return(sort(unique(microRNAome$cell_tissue)))
  if(!(x %in% microRNAome$cell_tissue)) return(NULL)
  x <- assay(microRNAome)[,microRNAome$cell_tissue==x,drop=FALSE]
  row.names(x) <- gsub("/.+","",(row.names(x)))
  x <- rowSums(x)
  x <- sort(log1p(10^6 * x/sum(x)),decreasing=TRUE)
  x[!duplicated(names(x))]
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
  x <- sort(miRNAexpMouse[,x], decreasing=TRUE)
  x[!duplicated(names(x))]
}

.addBestType <- function(x){
  stopifnot(!is.null(dim(x)))
  if(!is.null(x$best_stype)) return(x)
  if(length(w <- grep("[6-8]mer", colnames(x)))<2) return(x)
  x$best_stype <- apply(as.matrix(x[,w]), 1, FUN=function(y){
    if(length(i <- which(as.numeric(y)>0))==0) return(0)
    head(i,1)
  })
  x$best_stype <- factor(x$best_stype, levels=c(0L, rev(seq_along(w))), 
                         labels=c("no site",rev(gsub("Sites_","",colnames(x)[w]))))
  levels(x$best_stype) <- gsub("_","-",levels(x$best_stype))
  x
}

#' @importFrom data.table as.data.table setkey
.aggregateByFamilies <- function(ts, propTolerated=0.05){
  md <- metadata(ts)
  fam <- md$families
  if(is.null(md$families)) return(ts)
  prop <- sum(levels(ts$set) %in% md$families)
  if(prop>=(1-propTolerated)) return(ts)
  if(length(missng <- setdiff(levels(ts$set), names(md$families)))>0){
    warning("Some sets had no indicated family!")
    md$families <- c(
      setNames(as.character(md$families),names(md$families)),
      setNames(missng,missng))
  }
  levels(ts$set) <- as.character(md$families[levels(ts$set)])
  dt <- as.data.table(as.data.frame(ts))
  rm(ts)
  dt <- setkey(dt, set, feature)
  dt <- unique(dt[,c("set","feature","sites")], by=c("set","feature"))
  ts <- DataFrame(dt)
  rm(dt)
  metadata(ts) <- md
  ts
}



#' matchMirExpr
#' 
#' Tries to match input miRNA expression with given names, eventually adding the
#' -3p/-5p.
#'
#' @param me A named vector of miRNA expression, or data.frame with an 
#' 'expression' column
#' @param setsNames An optional vector of known sets names, or an actual target
#' annotation.
#' @param restrict Logical indicating whether to restrict `me` to elements 
#' included in `setNames` (after re-formating)
#'
#' @return A data.frame of miRNA expressions.
matchMirExpr <- function(me, setsNames=NULL, restrict=TRUE){
  if(is.vector(me)) me <- data.frame(row.names=names(me), expression=me)
  row.names(me) <- gsub("mir","miR", row.names(me))
  if(!any(grepl("miR-",row.names(me)))) return(me)
  w <- grep("-5p$|-3p$", row.names(me), invert=TRUE)
  if(length(w)==0) return(me)
  me2 <- data.frame(row.names=paste0(rep(row.names(me)[w],each=2),
                                     c("-5p","-3p")),
                    expression=rep(me$expression[w],each=2))
  me <- rbind(me,me2[setdiff(row.names(me2),row.names(me)),,drop=FALSE])
  if(!is.null(setsNames)){
    if(!is.vector(setsNames)){
      if(is(setsNames, "DFrame")){
        fams <- metadata(setsNames)$families
        setsNames <- as.character(unique(setsNames$set))
        if(!is.null(fams)) setsNames <- c(setsNames, unlist(fams), names(fams))
      }else if(is.list(setsNames)){
        setsNames <- names(setsNames)
      }else if(is.matrix(setsNames)){
        setsNames <- colnames(setsNames)
      }
    }
    nn <- setNames(gsub("-[1-9]-","-",row.names(me)), row.names(me))
    nn <- nn[nn %in% setsNames & !(nn %in% row.names(me)) & !duplicated(nn)]
    if(length(nn)>0){
      me2 <- data.frame(row.names=nn, expression=me[names(nn),"expression"])
      me <- rbind(me,me2)
    }
    if(restrict) me <- me[row.names(me) %in% setsNames,,drop=FALSE]
  }
  me
}

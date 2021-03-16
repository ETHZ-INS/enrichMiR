#' CDplot
#'
#' Creates a cumulative distribution plot, eventually handling division of 
#' values into bins according to another variable.
#'
#' @param ll A named list of numeric values, or a vector of values to be divided
#' @param by If `ll` is a list, `by` is ignored. If `ll` is a vector of values, 
#' `by` should be a numeric/logical/factor vector of the same length according 
#' to which the values should be divided.
#' @param k The number of divisions (ignored if `ll` is a list, or if `by` is
#' logical or a factor).
#' @param breaks The breaks to use (ignored if `ll` is a list). If NULL, will
#' be calculated.
#' @param sameFreq Logical; whether to calculate breaks so thate the groups have
#' the same frequency (rather than same width)
#' @param addN Logical; whether to add group sizes to their names
#' @param dig.lab Number of digits for automatically-generated breaks.
#' @param minN The minimum number of items per group (groups below this will be
#' merged if the breaks are numeric).
#' @param ... Passed to `geom_line` (can for instance be used for `size`, etc.)
#'
#' @return A ggplot.
#' @import ggplot2
#' @importFrom dplyr bind_rows
#' 
#' @export
CDplot <- function(ll, by=NULL, k=3, breaks=NULL, sameFreq=FALSE, addN=FALSE, 
                   dig.lab=NULL, minN=10, ...){
  library(ggplot2)
  if(!is.list(ll)){
    if(is.null(by)) stop("If `ll` is not already a list, `by` should be given.")
    if(length(by)!=length(ll)) stop("Lengths of ll and by differ.")
    w <- which(!is.na(by) & !is.na(ll))
    by <- by[w]
    ll <- ll[w]
    if(is.factor(by) || is.logical(by) || length(unique(by))<7){
      ll <- split(ll, by)
    }else{
      if(is.null(dig.lab)) dig.lab <- max(c(2,3-ceiling(log10(abs(mean(by))))))
      if(is.null(breaks)) breaks <- k
      if(sameFreq){
        k <- k+1
        breaks <- unique(quantile(by, prob=seq(from=0, to=1, length.out=k),na.rm=TRUE))
        if(length(breaks)<k)
          breaks <- unique(quantile(c(0,by[by!=0]), prob=seq(from=0, to=1, length.out=k),
                                    na.rm=TRUE))
        if(length(breaks)<k){
          desiredK <- k
          breaks <- 1
          while(length(breaks)<desiredK && k<100){
            k <- k
            breaks <- unique(quantile(by, prob=seq(from=0, to=1, length.out=k),
                                  na.rm=TRUE))
          }
        }
      }
      ll <- split(ll, cut(by, breaks, dig.lab=dig.lab))
    }
  }
  ll <- .mergeSmallerGroups(ll,minN=minN)
  d <- dplyr::bind_rows(lapply(ll, FUN=function(x){
    data.frame( y=(seq_along(x)-1)/(length(x)-1), x=sort(x) )
  }), .id="Sets")
  d$Sets <- factor(d$Sets, levels=unique(d$Sets))
  if(addN) levels(d$Sets) <- paste0(levels(d$Sets), " (n=",
                                    as.numeric(table(d$Sets)), ")")
  p <- ggplot(d, aes(x,y,colour=Sets)) + 
    geom_vline(xintercept=0, linetype="dashed") + geom_line(...) 
  p + ylab("Cumulative proportion")
}

#' CDplot2
#' 
#' A wrapper around \code{\link{CDplot}} to work directly with DEA results and
#' a target sets as inputs.
#'
#' @param dea A named numeric vector containing the signal for each feature, or
#' a DEA results data.frame (such as produced by major DEA packages, i.e. edgeR,
#' DESeq2 or limma).
#' @param sets The target sets object.
#' @param setName The name of the set to plot in `sets`.
#' @param k The number of groups to form.
#' @param by The variable by which to form the groups; either "sites" or "score"
#' @param sameFreq Logical; whether to divide so as to obtained groups with 
#' similar frequencies (rather than groups of similar width)
#' @param line.size Size of the line.
#' @param point.size Size of the points.
#' @param ... Any other argument passed to \code{\link{CDplot}}.
#'
#' @return A ggplot.
#' @export
CDplot2 <- function(dea, sets, setName, k=3, by=c("sites","score"), 
                    sameFreq=NULL, line.size=1.2, point.size=0.8, checkSynonyms=TRUE, ...){
  message("Preparing inputs...")
  dea <- .homogenizeDEA(dea)
  if(checkSynonyms) dea <- .applySynonyms(dea, sets)
  sets <- .list2DF(sets)
  sets <- sets[sets$set==setName,]
  if(is.null(sets$sites)) sets$sites <- 1
  by <- match.arg(by)
  if(!(by %in% colnames(sets))) stop("`by` not available in `sets`.")
  if(nrow(sets)==0) stop("setName not found in `sets`.")
  if(is.null(dim(dea))){
    if(!is.numeric(dea) || is.null(names(dea)))
      stop("`dea` should be a named numeric vector or a data.frame containing ",
           "the results of a differential expression analysis.")
  }else{
    dea <- .dea2sig(.homogenizeDEA(dea), "logFC")
  }
  if(!any(sets$feature %in% names(dea)))
    stop("There seems to be no overlap between the rows of `dea` and the ",
         "features of `sets`.")
  if(is.null(sameFreq)) sameFreq <- by=="score"
  
  by <- rowsum(sets[[by]], sets$feature)[,1]
  by2 <- by[names(dea)]
  by2[is.na(by2)] <- 0
  if(k==2){
    ll <- split(dea, by2!=0)
    names(ll) <- c("non-targets","targets")
    p <- CDplot(ll, size=line.size, ...)
  }else{
    p <- CDplot(dea, by=by2, k=k, sameFreq=sameFreq, size=line.size, ...)
  }
  if(point.size>0) p <- p + geom_point(size=point.size)
  p + xlab("logFC") + ggtitle(setName)
}

.mergeSmallerGroups <- function(ll, minN=10){
  if(all(grepl("^[\\(\\[\\)\\]].+,.+[\\(\\[\\)\\]]$", names(ll), perl=TRUE))){
    lmf <- function(x){
      x <- strsplit(x,",")
      paste(x[[1]][[1]],x[[2]][[2]],sep=",")
    }
  }else if(all(grepl("^[0-9]+$", names(ll)))){
    lmf <- function(x) paste(x,sep="-")
  }else{
    return(ll)
  }
  while(length(ll)>2 && any(lengths(ll)<minN)){
    w <- which.min(lengths(ll))
    if(w==1){ w2 <- 2 }
    else if(w==length(ll)){ w2 <- length(ll)-1 }
    else{
      w2 <- w+which.min(lengths(ll)[w-1],lengths(ll)[w+1])*2-3
    }
    ll[[w]] <- sort(c(ll[[w]],ll[[w2]]))
    names(ll[[w]]) <- lmf(names(ll)[sort(c(w,w2))])
    ll[[w2]] <- NULL
  }
  ll   
}





#' enrichPlot
#'
#' Plots the results of an miRNA enrichment analysis
#'
#' @param res The results of an enrichment analysis, as obtained by the 
#' `getResults` function.
#' @param enr.field The column determining the x axis.
#' @param size.field The column determining the size of the circles.
#' @param col.field The column determining the color of the circles.
#' @param sig.field The column to use as significance or y axis (default 'FDR').
#' @param max.sig.thres The maximum FDR below which to plot. Alternatively, if
#' a number > 1, indicates the number of top points to plot.
#' @param min.enr.thres The minumum (absolute) fold-enrichment above which to plot.
#' @param label.sig.thres The FDR threshold above which to plot the names of the
#'  enriched sets (default 0.05).
#' @param label.enr.thres The (absolute) enrichment threshold above which to 
#' plot the names of the enriched sets.
#' @param label.field The field to use as labels (defaults to row.names)
#' @param maxLabels The maximum number of labels to plot (default 10).
#' @param opacity The opacity of the bubbles, ranging from 0 to 1 (default 0.5).
#' @param repel Logical; whether to plot labels with `ggrepel`. This is the 
#' default behaviour, but should be turned off if you plan to pass the result
#' to `ggplotly`. 
#'
#' @export
#' @import ggplot2
enrichPlot <- function( res,
                        enr.field=c("enrichment","normalizedEnrichment","beta","coefficient"),
                        size.field=c("overlap", "set_size"),
                        col.field=NULL,
                        sig.field="FDR", 
                        max.sig.thres=100, 
                        min.enr.thres=NULL, 
                        label.sig.thres=0.05, 
                        label.enr.thres=0, 
                        label.field=NULL, 
                        maxLabels=10, 
                        opacity=0.5,  
                        repel=TRUE ){
  if(length(enr.field <- head(intersect(enr.field,colnames(res)),n=1))==0)
    stop("`enr.field` not found.")
  if(length(sig.field <- head(intersect(sig.field,colnames(res)),n=1))==0)
    stop("`sig.field` not found.")
  if(!is.null(size.field) && 
     length(size.field <- head(intersect(size.field,colnames(res)),n=1))==0)
    stop("`size.field` not found.")
  if(!is.null(col.field) && 
     length(col.field <- head(intersect(col.field,colnames(res)),n=1))==0)
    stop("`col.field` not found.")
  if(is.null(label.field)){
    res$set <- row.names(res)
    label.field <- "set"
  }else{
    if(length(label.field<- head(intersect(label.field,colnames(res)),n=1))==0)
      stop("`label.field` not found.")
  }
  if(!is.null(res[["overlap"]])) res[["overlap"]] <- as.numeric(res[["overlap"]])
  if(!is.null(res[["set_size"]])) res[["set_size"]] <- as.numeric(res[["set_size"]])
  
  if(!is.null(min.enr.thres)) res <- res[res[[enr.field]]>=min.enr.thres,]
  res <- res[order(res[[sig.field]]),]
  if(max.sig.thres>1){
    res <- res[head(order(res[[sig.field]]),n=max.sig.thres),]
  }else{
    res <- res[which(res[[sig.field]] <= max.sig.thres),]
  }
  if(nrow(res)<=2) stop("Not enough data to plot using the given thresholds...")
  w <- which(abs(res[[enr.field]])>=label.enr.thres & 
               res[[sig.field]]<label.sig.thres)
  w <- head(w,n=maxLabels)
  sig.field2 <- paste0("-log10(",sig.field,")")
  res[[sig.field2]] <- -log10(res[[sig.field]])
  ll <- list(label=label.field, x=enr.field, y=sig.field2)
  if(!is.null(size.field)) ll$size <- size.field
  if(!is.null(col.field)) ll$colour <- col.field
  for(f in setdiff(colnames(res), unlist(ll))) ll[[f]] <- f
  p <- ggplot(res, do.call(aes_string, ll)) + geom_point(alpha=opacity)
  if(repel){
    p <- p + geom_text_repel(data=res[w,])
  }else{
    p <- p + geom_text(data=res[w,])
  }
  p
}


#' breakStrings
#'
#' Breaks each (long) component of a character vector into two lines.
#'
#' @param x A character vector.
#' @param minSizeForBreak The minimum number of characters for a string to be 
#' broken onto two lines (default 20).
#' @param lb The linebreak character to use (default \\n).
#'
#' @export
breakStrings <- function (x, minSizeForBreak = 20, lb = "\n") {
    sapply(x, minSizeForBreak = minSizeForBreak, lb = lb, FUN = function(x, 
        minSizeForBreak, lb) {
        if (nchar(x) <= minSizeForBreak) 
            return(x)
        g <- gregexpr(" ", x)[[1]]
        if (length(g) == 0) 
            return(x)
        if (length(g) == 1 & all(g == -1)) 
            return(x)
        mid <- nchar(x)/2
        mid <- g[order(abs(g - mid))[1]]
        substr(x, mid, mid) <- lb
        return(x)
    })
}

#' CDplot
#'
#' Creates a cumulative distribution plot, eventually handling division of 
#' values into bins according to another variable.
#'
#' @param x A named list of numeric values, or a vector of values to be divided
#' @param by If `ll` is a list, `by` is ignored. If ll vector of values, `by` 
#' should be a vector of the same length according to which the values should
#' be divided.
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
      if(sameFreq)
        breaks <- unique(quantile(by, prob=seq(from=0, to=1, length.out=k+1),
                                  na.rm=TRUE))
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
#' @param res The results of an enrichment analysis, as obtained by the `enrichMiR.results` function.
#' @param size.field The column determining the size of the bubbles (default "overlap").
#' @param col.field The column determining the color of the bubbles (default "expression" if present, otherwise "overlap").
#' @param sig.field The column to use as significance (default 'FDR').
#' @param min.sig.thres The minimum FDR above which to plot (default 0.2).
#' @param min.enr.thres The minumum fold-enrichment above which to plot (default 1).
#' @param label.sig.thres The FDR threshold above which to plot the names of the enriched families (default 0.05).
#' @param label.enr.thres The enrichment threshold above which to plot the names of the enriched families  (default 2).
#' @param label.field The field to use as labels (default 'family')
#' @param maxLabels The maximum number of labels to plot (default 10).
#' @param opacity The opacity of the bubbles, ranging from 0 to 1 (default 0.5).
#' @param main The title of the plot.
#'
#' @export
enrichPlot <- function( res, 
                        size.field="overlap", 
                        col.field=NULL,
                        sig.field="FDR", 
                        min.sig.thres=0.2, 
                        min.enr.thres=1, 
                        label.sig.thres=0.05, 
                        label.enr.thres=2, 
                        label.field="family", 
                        maxLabels=10, 
                        opacity=0.5,  
                        main=""){
    if(is(res,"enrichMiR")){
      res <- enrichMiR.results(res)
      if(!("enrichment" %in% colnames(res))){
        if("EN.combined.enrichment" %in% colnames(res)){
          res$enrichment <- res$EN.combined.enrichment
          res$overlap <- res$EN.combined.overlap
          min.enr.thres <- -Inf
          xlab <- "Combined log2(fold-enrichment)"
        }else{
          res$enrichment <- apply(abs(log2(as.matrix(res[,grep("enrichment"),drop=F])+0.05)),1,FUN=median)
          xlab <- "Median absolute log2(fold-enrichment)"
        }
      }
      res$medianP <- apply(res[,grep("pvalue",colnames(res)),drop=F],1,FUN=median)
      sig.field <- "medianP"
    }else{
      xlab <- "Fold enrichment"
      if(!("enrichment" %in% colnames(res))){
        if("normalizedEnrichment" %in% colnames(res)){
          res$enrichment <- res$normalizedEnrichment
          xlab <- "Normalized enrichment"
        }else{
          if("logFC" %in% colnames(res)){
            res$enrichment <- res$logFC
            xlab <- "logFC per unit"
          }else{
            stop("These results do not seem to contain an enrichment value.")
          }
        }
      }
    }
    if(is.null(col.field)){
      if("expression" %in% colnames(res)){
        col.field <- "expression"
      }else{
        if("overlap" %in% colnames(res)){
          col.field <- "overlap"
        }else{
          col.field <- "enrichment"
        }
      }
    }
    suppressPackageStartupMessages(library(plotly))
    if("features" %in% colnames(res)){
      if(is(res$features, "CharacterList")){
        res <- res[which(!duplicated(paste(res$features, collapse=""))),]
      }else{
        res <- res[which(!duplicated(res$features)),,drop=F]
      }
    }
    if(nrow(res)<=2) stop("Not enough data to plot! Perhaps check the thresholds being used...")
    res <- res[which( res$enrichment>=min.enr.thres & res[[sig.field]] <= min.sig.thres),]
    w <- which(res$enrichment>label.enr.thres & res[[sig.field]]<label.sig.thres)

    if("expression" %in% colnames(res)){
      lab <- paste0(row.names(res)," (expr~",round(res$expression,1),")")
    }else{
      lab <- row.names(res)
    }
    lab <- paste0(lab,"\n",res$miRNAs,"\n",
                  round(res$enrichment,2),ifelse(xlab=="Fold enrichment","-"," "),xlab,"\n",
                  sig.field,"~",format(res[[sig.field]],digits=2))
    if("overlap" %in% colnames(res)) lab <- paste0(lab, "\n", res$overlap, " target features")
    if("BS.in" %in% colnames(res)) lab <- paste0(lab, "\n", res$BS.in, " binding sites")
    
    if(length(w)>0){
        if(length(w)>maxLabels) w <- w[seq_len(maxLabels)]
        if(label.field=="family" | label.field=="row.names"){
          lab2 <- row.names(res)[w]
        }else{
          lab2 <- res[[label.field]][w]
        }
        a <- list( x = res$enrichment[w], y=-log10(res[[sig.field]][w]), 
                   text=breakStrings(lab2),
                   showarrow = FALSE, ax = 0, ay = 0)
    }

    ml <- list(opacity=opacity)
    if(!is.null(size.field)) ml$size <- as.formula(paste0("~",size.field))
    p <- plot_ly( res, 
                  x=res$enrichment, 
                  y=-log10(res[[sig.field]]),
                  type='scatter',
                  mode='markers', 
                  text=lab, 
                  hoverinfo='text', 
                  color=as.formula(paste0("~",col.field)), 
                  marker=ml ) %>% 
        layout( title = main,
                xaxis = list(showgrid = FALSE, title=xlab),
                yaxis = list(showgrid = FALSE, title=paste0("-log10(",sig.field,")"))  )

    if(length(w)>0) p <- p %>% layout(annotations=a, font=list(color="black",size=13))

    p
}


#' breakStrings
#'
#' Breaks each (long) component of a character vector into two lines.
#'
#' @param x A character vector.
#' @param minSizeForBreak The minimum number of characters for a string to be broken onto two lines (default 20).
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

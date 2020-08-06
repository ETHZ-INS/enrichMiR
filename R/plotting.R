#' CDplot
#'
#' Plots the cumulative foldchange distribution for a given miRNA family or lists of foldchanges.
#'
#' @param x Either and object of class `enrichMiR`, or a list containing vectors of foldchanges.
#' @param family If `x` is an object of class `enrichMiR`, the family to be plotted.
#' @param col The vector of colors to be used for the different distributions (defaults to R default colors).
#' @param xlim The x axis range to be displayed, default `c(-2,2)`.
#' @param xlab The x axis label, defaults to "log2(foldchange)".
#' @param main The plot's title, defaults to the family (if given).
#' @param sub The plot's subtitle, defaults to the representative miRNA if the `x` is an object of calss `enrichMiR`.
#' @param splitBy Split the genes by number of binding sites (`splitBy='sites'`, default), by targetscan score 
#' (`splitBy='score'`), or simply into target/non-target (`splitBy=NULL`).
#' @param ... Any other arguments passed to the plot function.
#'
#' @export
CDplot <- function(x, family=NULL, cols=1:8, xlim=c(-2,2), xlab="log2(foldchange)", main="", sub=NULL, splitBy="sites", ...){
    if(is(x,"enrichMiR")){
        if(is.null(sub) && !is.na(x@families)) sub <- toString(names(x@families)[which(x@families==family)])
        if(!any(x@TS$family==family)) stop("miRNA family not found!")
        if(is.null(main) || main=="") main <- family
        m <- merge(x@TS[which(x@TS$family==family),], x@DEA, by.x="feature", by.y="row.names",all.y=T)
        m$sites[which(is.na(m$sites))] <- 0
        m$score[which(is.na(m$score))] <- 0
        if(is.null(splitBy)){
          w <- which(m$sites>0)
          ll <- list('non-targets'=m$feature[-w], targets=m$feature[w])
        }else{
          if(splitBy=="score"){
            w <- which(m$score>0)
            q <- unique(quantile(m$score[w],probs=c(0,0.333,0.666,1)))
            print(q)
            ll <- split(m$feature[w], cut(m$score[w],q,include.lowest = T))
            print(names(ll))
            ll <- c(list('non-targets'=m$feature[-w]),ll)
          }else{
            tt <- table(m$sites)
            bk1 <- unique(c(0,as.numeric(names(tt)[which(tt>8)]),Inf))
            bk <- cut(m$sites, breaks=bk1, include.lowest=T, right=F)
            bk1[length(bk1)-1] <- paste0(">",bk1[length(bk1)-1])
            levels(bk) <- paste(bk1[-length(bk1)], c("site",rep("sites",length(bk1)-2)))
            ll <- split(m$feature, bk)
          }
        }
        ll <- lapply(ll, dea=x@DEA, FUN=function(x,dea){ dea[x,"logFC"] })
        CDplot(ll, main=main, sub=sub, xlim=xlim, xlab=xlab, cols=cols, ...)
    }else{
        x2 <- lapply(x,FUN=ecdf)
        p <- sapply(x[-1],a=x[[1]],FUN=function(x,a){ suppressWarnings(ks.test(x,a)$p.value) })
        plot(x2[[1]],col=cols[1],xlab=xlab,ylab="Cumulative proportion",xlim=xlim,main=main,sub=sub,...)
        abline(v=0,lty="dashed",col="grey")
        for(i in 2:length(x)) lines(x2[[i]],col=cols[i])
        legend("bottomright",bty="n",lwd=3,col=cols,legend=paste(names(x), c("",paste0(" (p~",sapply(p,digits=2,FUN=format),")"))))
    }
}


CDplot <- function(ll, by=NULL, k=5, ...){
  if(!is.list(ll)){
    if(is.null(by)) stop("If `ll` is not already a list, `by` should be given.")
    if(length(by)!=length(ll)) stop("Lengths of ll and by differ.")
    ll <- split(ll, cut(by, k, dig.lab=3-ceiling(log10(abs(mean(by))))))
  }
  p <- format(suppressWarnings(ks.test(ll[[1]], rev(ll)[[1]])$p.value), digits=2)
  message("KS p-value between first and last sets:\n", p)
  d <- dplyr::bind_rows(lapply(ll, FUN=function(x){
    data.frame( y=seq_along(x)/length(x),
                x=sort(x) )
  }), .id="Genesets")
  d$Genesets <- factor(d$Genesets, levels=unique(d$Genesets))
  p <- ggplot(d, aes(x,y,colour=Genesets)) + 
    geom_vline(xintercept=0, linetype="dashed") + geom_line(...)
  p + ylab("Cumulative proportion")
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

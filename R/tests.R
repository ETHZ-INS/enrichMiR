#' overlap
#'
#' Traditional overlap (Fisher test) between differentially-expressed features and miRNA targets.
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param sets A data.frame of annotation, with at least the following columns: `set`, `feature`.
#'
#' @return a data.frame.
#'
#' @export
overlap <- function(signal, sets){
  tested <- names(signal)
  significant <- names(signal)[signal]
  significant <- intersect(significant,tested)
  res <- vapply( split(sets$feature,sets$set), set1=significant, 
                FUN.VALUE=numeric(5), FUN=function(set2,set1){
    expected <- length(set1)*length(set2)/length(tested)
    ov <- length(intersect(set1,set2))
    c(annotated=length(set2),
      overlap=ov,
      enrichment=round(log2((ov+0.5)/expected),2),
      under.pvalue=.overlap.prob(set1,set2,tested,lower=T),
      over.pvalue=.overlap.prob(set1,set2,tested)
    )
  })
  res <- as.data.frame(t(res))
  res$FDR <- p.adjust(res$over.pvalue,method="fdr")
  res[order(res$FDR,res$over.pvalue),]
}


#' plMod
#'
#' @param signature A named numeric vector
#' @param sets A data.frame with at least the following columns: 'set', 
#' 'feature', and the variable specified by `var`.
#' @param var The independent variable to use, either 'sites' or 'score'.
#' @param correctForLength Logical; whether to correct for length (i.e. total
#'  sites). Defaults to TRUE if `var='sites'`, FALSE otherwise.
#'
#' @return a data.frame
#' @export
plMod <- function(signature, sets, var="sites", correctForLength=NULL){
  if(is.null(correctForLength)) 
    correctForLength <- (var=="sites" && 'sites' %in% colnames(sets))
  if(correctForLength){
    ag <- aggregate(sets$sites, by=list(feature=sets$feature), FUN=sum)
    cfl <- ag[,2]
    names(cfl) <- ag[,1]
    cfl <- cfl[names(signature)]
    cfl[which(is.na(cfl))] <- 0
    names(cfl) <- names(signature)
  }else{
    cfl <- NULL
  }
  sets$family <- as.character(sets$set)
  res <- as.data.frame(t(vapply(split(sets,sets$set), fcs=signature, cfl=cfl, 
                  FUN.VALUE=numeric(2), FUN=function(x,fcs, minSize, cfl){
    x <- x[!duplicated(x),]
    row.names(x) <- x$feature
    x2 <- x[names(fcs),var]
    x2[which(is.na(x2))] <- 0
    if(is.null(cfl)){
      mod <- try(lm(fcs~x2+0),silent=T)
    }else{
      mod <- try(lm(fcs~x2+cfl+0),silent=T)
    }
    if(!is(mod,"try-error"))
      return(c(coef(mod)["x2"], summary(aov(mod))[[1]]["x2","Pr(>F)"]))
    return(rep(NA_real_,2))
  })))
  colnames(res) <- c("coefficient","pvalue")
  res$FDR <- p.adjust(res$pvalue)
  res[order(res$FDR,res$pvalue),]
}

modsites <- function(x, sets, correctForLength=TRUE, ...){
  plMod(x, sets, correctForLength=correctForLength, ...)
}

modscore <- function(x, sets, ...) plMod(x, sets, var="score", ...)


#' woverlap
#'
#' Weighted hypergeometric test adjusted for total number of miRNA bindings sites. (This is done through the `goseq` package)
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param sets A data.frame with at least the following columns: 'set', 'feature'.
#' @param method Method for handling then length bias, default "Wallenius". See `?goseq` for more detail.
#'
#' @return a data.frame.
#'
#' @importFrom goseq nullp goseq
#' @export
woverlap <- function(signal, sets, method="Wallenius"){
  if(is.null(sets$sites)){
    bd <- as.numeric(table(sets$feature)[names(signal)])
    names(bd) <- names(signal)
    bd[is.na(bd)] <- 0
  }else{
    bd <- sapply(signal, FUN=function(x) 0)
    ag <- aggregate(sets$sites, by=list(gene=sets$feature), FUN=sum)
    bd[ag$gene] <- ag$x
  }
  np <- nullp(signal[names(bd)], bias.data = bd, plot.fit=FALSE)
  g2c <- split(sets$feature, sets$set)

  res <- goseq(np, gene2cat=g2c, method=method, use_genes_without_cat=TRUE)
  row.names(res) <- res[,1]
  res <- res[,-1]
  colnames(res) <- c("over.pvalue","under.pvalue","overlap","numInCat")

  # TEMPORARY enrichment value
  significant <- names(signal)[signal]
  res$enrichment <- round(log2(res$overlap/(length(significant)*(res$numInCat/length(signal)))),2)
  
  # temporary fix of pvalues >1 or <0
  res$over.pvalue[res$over.pvalue > 1] <- 1
  res$under.pvalue[res$under.pvalue > 1] <- 1
  res$over.pvalue[res$over.pvalue < 0] <- 0
  res$under.pvalue[res$under.pvalue < 0] <- 0
  
  res$FDR <- p.adjust(res$over.pvalue, method="fdr")
  colnames(res)[4] <- "annotated"
  res
}

#' siteoverlap
#'
#' Applies Fisher's test to the number of miRNA binding sites among a set of features (vs all other binding sites in that set of features).
#' Of note, this application violates the assumptions of Fisher's exact test.
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param sets A data.frame with at least the following columns: 
#' 'set', 'feature', 'sites'.
#'
#' @return a data.frame.
#'
#' @export
siteoverlap <- function(signal, sets){
  tested <- names(signal)
  significant <- names(signal)[signal]
  allBS.bg <- sum(sets[which(!(sets$feature %in% significant)),"sites"])
  allBS.sig <- sum(sets[which(sets$feature %in% significant),"sites"])
  res <- t(vapply(split(sets,sets$set), set1=significant, bs.sig=allBS.sig, 
                  bs.bg=allBS.bg, FUN.VALUE=numeric(9), 
                  FUN=function(x,set1,bs.sig,bs.bg){
    x <- as.data.frame(x)
    w <- which(as.character(x$feature) %in% set1)
    xin <- sum(x[w,"sites"])
    xout <- sum(x[which(!(as.character(x$feature) %in% set1)),"sites"])
    mm <- matrix(round(c(xin,bs.sig-xin,xout,bs.bg-xout)),nrow=2)
    p1 <- fisher.test(mm,alternative="greater")$p.value[[1]]
    p2 <- fisher.test(mm,alternative="less")$p.value[[1]]
    c(  annotated=nrow(x),
        overlap=length(unique(x[w,"feature"])),
        BS.in=xin,
        otherBS.in=bs.sig-xin,
        BS.inBG=xout,
        enrichment=log2((xin/(bs.sig-xin))/(xout/(bs.bg-xout))),
        otherBS.inBG=bs.bg-xout,
        under.pvalue=p2,
        fisher.pvalue=p1
    )
  }))
  res <- as.data.frame(res)
  res <- res[order(res$fisher.pvalue),]
  for(i in 1:5) res[[i]] <- as.integer(res[[i]])
  res$FDR <- p.adjust(res$fisher.pvalue,method="fdr")
}

#' ks
#'
#' enrichment analysis using a Kolmogorov-Smirnov test on the signal.
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param sets A data.frame with at least the following columns: 'set', 'feature'.
#'
#' @return a data.frame.
ks <- function(signal, sets){
  res <- t(vapply(split(sets$feature,sets$set), FUN.VALUE=numeric(2), FUN=function(x){
    ks <- suppressWarnings(try(ks.test(signal[x], signal[setdiff(names(signal),x)])$p.value, silent=T))
    if(is(ks,"try-error")) ks <- NA
    c( annotated=length(set),
       ks.pvalue=ks )
  }))
  res <- as.data.frame(res)
  res$FDR <- p.adjust(res$ks.pvalue,method="fdr")
  res[order(res$FDR,res$ks.pvalue),]
}

#' mw
#'
#' miRNA targets enrichment analysis using a Mann-Whitney / Wilcoxon test on the targets' foldchanges.
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param sets A data.frame with at least the following columns: 'set', 'feature'.
#'
#' @return a data.frame.
mw <- function(signal, sets){
  signal <- names(signal)[signal]
  res <- t(vapply(split(sets$feature,sets$set), FUN.VALUE=numeric(2), FUN=function(x){
    set <- intersect(unique(as.character(x$feature)),names(signal))
    mw <- suppressWarnings(try(wilcox.test(signal[set], signal[setdiff(names(signal),set)])$p.value, silent=T))
    if(is(mw,"try-error")) mw <- NA
    c( annotated=length(set), ks.pvalue=ks )
  }))
  res <- as.data.frame(res)
  res$FDR <- p.adjust(res$mw.pvalue,method="fdr")
  res[order(res$FDR,res$mw.pvalue),]
}


#' gsea
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param sets A data.frame with at least the following columns: 'set', 'feature'.
#' @param maxSize The maximum number of elements in a set to be tested (default 500). If the test takes too long to run, consider setting this.
#' @param nperm The number of permutations, default 2000. The more permutations, the more precise the p-values, in particular for small p-values.
#'
#' @return a data.frame.
#'
#' @importFrom fgsea fgsea
#' @export
gsea <- function(signal, sets, maxSize=300, nperm=2000, ...){
  sets <- lapply(split(sets$feature,sets$set), tested=names(signal), 
                 FUN=function(x,tested){ intersect(unique(x),tested) })
  res <- fgsea(sets, signal, nperm, minSize=4, maxSize=maxSize, ...)
  res <- res[order(res$padj,res$pval),]
  colnames(res)[1:5] <- c("family","pvalue","FDR","ES","normalizedEnrichment")
  colnames(res)[8] <- "features"
  row.names(res) <- res[,1]
  return(res[,-1])
}


#' geneBasedTest
#'
#' Applies Fisher's test to the number of miRNA binding sites of each differentially-expressed gene
#'
#' @param features The set of differentially-expressed features, or features of interest.
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#'
#' @return a named vector of p-values.
#'
#' @importFrom aggregation fisher
#' @export
geneBasedTest <- function(features, TS){
  TS$family <- as.character(TS$family)
  features <- as.character(features)[which(features %in% TS$feature)]
  tmp <- aggregate(TS$sites,by=list(family=TS$family),FUN=sum)
  bgs <- tmp[,2]
  names(bgs) <- tmp[,1]
  TS <- TS[which(TS$feature %in% features),]
  
  res <- sapply(split(TS[,c("family","sites")],as.character(TS$feature),drop=F), fams=unique(TS$family), bgs=bgs, FUN=function(x, fams, bgs){
    p <- rep(1,length(fams))
    names(p) <- fams
    tbs.in <- sum(x$sites)
    tbs <- sum(bgs)
    x <- as.data.frame(x)
    for(f in 1:nrow(x)){
      y <- x[f,"sites"]
      b <- bgs[x[f,"family"]]
      mm <- matrix(c(y, sum(x$sites)-y, b, sum(bgs)-b),nrow=2)
      r <- try(fisher.test(mm, alternative="greater")$p.value,silent=T)
      p[x[f,"family"]] <- ifelse(is(r,"try-error"),NA,r)
    }
    p
  })
  res <- apply(res,1,FUN=function(x){ x <- x[!is.na(x)]; if(length(x)==0) return(NA); fisher(x)})
  return(res)
}


.overlap.prob <- function (set1, set2, universe, lower = F){
  set1 <- as.character(set1)
  set2 <- as.character(set2)
  if (class(universe) == "character") {
    set1 <- intersect(set1, universe)
    set2 <- intersect(set2, universe)
    universe <- length(unique(universe))
  }
  set1 <- unique(set1)
  set2 <- unique(set2)
  ov <- sum(set1 %in% set2)
  phyper(max(0, ov - 1), length(set1), universe - length(set1), 
         length(set2), lower.tail = lower)
}



.censorScore <- function(x){
  x <- 0.1-x
  if(length(w <- which(x>0.9))>0)
    x[w] <- 0.9+ecdf(x[w])(x[w])/10
  x
}

TS2regulon <- function(x, likelihood="score"){
  x <- as.data.frame(x)
  if(likelihood=="score"){
    x$likelihood <- .censorScore(x[[likelihood]])
  }else{
    x$likelihood <- x[[likelihood]]
  }
  lapply(split(x,x$set, drop=TRUE), FUN=function(x){
    y <- list(  tfmode=rep(-1,nrow(x)),
                likelihood=x$likelihood )
    lapply(y, FUN=function(a){
      names(a) <- x$feature
      a
    })
  })
}

#' areamir
#'
#' analytic Rank-based Enrichment Analysis using a conversion of the scores as 
#' weights.
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param sets A data.frame with at least the following columns: 'set', 'feature', 'score'.
#'
#' @return a data.frame.
#'
#' @export
areamir <- function(signal, sets, ...){
  vi <- viper::msviper(signal, regulon=TS2regulon(as.data.frame(sets)), ..., verbose=FALSE)
  vi2 <- as.data.frame(vi$es[c("nes","size","p.value","nes.bt")])
  colnames(vi2)[2:3] <- c("annotated","pvalue")
  vi2$FDR <- p.adjust(vi2$pvalue)
  vi2 <- vi2[order(vi2$pvalue),]
  vi2
}

#' regmir
#'
#' miRNA enrichment analysis using regularized regression
#'
#' @param signal A vector of logical or numeric values, with gene symbols as names.
#' @param sets A data.frame with at least the following columns: 
#' 'set', 'feature', and (if binary=FALSE) 'score'.
#' @param binary Logical; whether to consider target prediction as binary.
#' @param alpha elastic net mixing param (0=ridge, 1=lasso)
#' @param do.plot Logical, whether to plot coefficients against lambda
#' @param use.intercept Logical, whether to use an intercept in the model.
#' @param keepAll Logical, whether to return all families.
#'
#' @return A DataFrame.
#' @import glmnet zetadiv S4Vectors 
#' @export
regmir <- function(signal, sets, binary=NULL, alpha=1, do.plot=FALSE, use.intercept=FALSE,
                   keepAll=TRUE){
  if(is.null(names(signal))) stop("`signal` should be a named vector!")
  suppressPackageStartupMessages(c(
    library(glmnet),
    library(zetadiv)
    ))
  if(is.null(binary)) binary <- "score" %in% colnames(sets)
  # prepare the target matrix
  if(binary){
    bm <- sapply(split(sets$feature, sets$set), FUN=function(x) names(signal) %in% x)
  }else{
    TS2 <- aggregate(sets[,"score",drop=FALSE], by=as.data.frame(sets[,c("set","feature")]), FUN=min)
    TS2 <- split(TS2[,c("score","feature")], TS2$set)
    bm <- sapply(TS2, FUN=function(x){
      row.names(x) <- x$feature
      -1*x[names(signal),"score"]
    })
    bm[is.na(bm)] <- 0
    colnames(bm) <- names(TS2)
  }
  bm <- bm[,colSums(bm)>0]
  
  # regularized regression with cross-validation
  if(isLogistic <- is.logical(signal)){
    fits <- cv.glmnet(bm, signal, standardize=FALSE, alpha=alpha, family="binomial", lower.limits=0 )
  }else{
    fits <- cv.glmnet(bm, signal, standardize=FALSE, alpha=alpha, family="gaussian")
  }
  
  if(do.plot){
    layout(matrix(1:2,nrow=1))
    plot(fits)
    plot(fits$glmnet.fit, label=TRUE)
  }
  
  # we extract the miRNAs selected by the best most regularized glmnet fit:
  co <- coef(fits, fits$lambda.1se)
  co <- row.names(co)[co[,1]!=0][-1]
  
  if( length(co)==0 && fits$lambda.min!=fits$lambda.1se ){
    # if no coefficient was selected, we use the minimum lambda
    co <- coef(fits, fits$lambda.min)
    co <- row.names(co)[co[,1]!=0][-1]
  }
  
  signald <- data.frame( y=signal, bm[,co,drop=FALSE] )
  
  # new fit to get significance estimates
  if(use.intercept){
    form <- y~.
  }else{
    form <- y~0+.
  }
  if(isLogistic){
    mod <- glm.cons( form, data=signald, family="binomial", cons=1, cons.inter=-1)
  }else{
    mod <- lm( form, data=signald )
  }
  
  # we extract the coefficients and p-values, and reorganize the output:
  res <- coef(summary(mod))
  res <- res[order(res[,4]),,drop=FALSE]
  colnames(res) <- c("beta","stderr",ifelse(isLogistic,"z","t"),"pvalue")
  res <- res[grep("^\\(Intercept\\)$|FALSE$", row.names(res), invert=TRUE),,drop=FALSE]
  row.names(res) <- gsub("TRUE","",row.names(res))
  
  if(keepAll){
    co2 <- sort(apply(fits$glmnet.fit$beta,1,FUN=function(x){
      if(!any(x!=0)) return(Inf)
      which(x!=0)[1]
    }))
    co2 <- co2[grep("^\\(Intercept\\)$|FALSE$", names(co2), invert=TRUE)]
    names(co2) <- gsub("TRUE","",names(co2))
    co2 <- co2[setdiff(names(co2),row.names(res))]
    co2 <- data.frame(row.names=names(co2), beta=rep(NA_real_,length(co2)),
                      stderr=NA_real_, z=NA_real_, pvalue=1)
    if(!isLogistic) colnames(co2)[3] <- "t"
    res <- rbind(res,co2)
  }
  
  # res <- DataFrame(res)
  # we adjust using all features as number of comparisons
  if(nrow(res)>0){
    res$FDR <- p.adjust(res$pvalue, n=ncol(bm))
    # res$features <- CharacterList(lapply(split(TS$feature, TS$family)[row.names(res)], 
    #                                 y=names(signal)[signal], FUN=intersect))
  }
  res
}

regmirb <- function(signal, ...){
  regmir(signal, binary=TRUE, ...)
}
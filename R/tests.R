#' overlap
#'
#' Traditional overlap (Fisher test) between a logical signal (e.g. 
#' differentially-expressed features) and the elements of sets.
#'
#' @param signal A named logical vector indicating which features are 
#' differentially-expressed
#' @param sets A data.frame of annotation, with at least the following columns: 
#' `set`, `feature`. Alternatively, a sparse logical matrix, with sets as 
#' columns and features as rows (dimensions must be named).
#' @param alternative Test alternative (defaults to 'greater' to test for
#' over-representation)
#'
#' @return a data.frame.
#'
#' @export
overlap <- function(signal, sets, alternative=c("greater","less","two.sided")){
  alternative <- match.arg(alternative)
  if(!.checkSets(sets, matrixAlternative="logical")){
    sets <- .aggregateByFamilies(sets)
    if("sites" %in% colnames(sets)){
      sets$sites <- sets$sites > 0
      sets <- .setsToScoreMatrix(signal, sets, column="sites", keepSparse=TRUE)
    }else{
      sets <- .setsToScoreMatrix(signal, sets, column=NULL, keepSparse=TRUE)
    }
  }else{
    sets <- sets[row.names(sets) %in% names(signal),]
  }
  signal <- signal[row.names(sets)]
  sigAndSet <- colSums(sets[which(signal),,drop=FALSE])
  notSigAndSet <- colSums(sets[which(!signal),,drop=FALSE])
  sigAndNotSet <- sum(signal)-sigAndSet
  notNot <- nrow(sets)-sigAndSet-notSigAndSet-sigAndNotSet
  expected <- sum(signal)*(sigAndSet+notSigAndSet)/nrow(sets)
  res <- data.frame( overlap=sigAndSet,
                     enrichment=round(log2((sigAndSet+0.25)/(expected+0.25)),3),
                     pvalue=fisher.test.p(sigAndSet,notSigAndSet,
                                               sigAndNotSet,notNot,
                                               alternative=alternative) )
  res$FDR <- p.adjust(res$pvalue, method="fdr")
  res[order(res$FDR,res$pvalue),]
}

#' siteoverlap.old
#'
#' Applies Fisher's test to the number of binding sites among a set of features 
#' (vs all other binding sites in that set of features).
#' This is an old apply-based implementation. For a much faster vectorial 
#' implementation, see `siteoverlap2`
#'
#' @param signal A named logical vector indicating which features are 
#' differentially-expressed
#' @param sets A data.frame with at least the following columns: 
#' 'set', 'feature', 'sites'.
#' @param alternative Test alternative (defaults to 'greater' to test for
#' over-representation)
#'
#' @return a data.frame.
#'
#' @export
siteoverlap.old <- function(signal, sets){
  sets <- .aggregateByFamilies(sets)
  tested <- names(signal)
  significant <- names(signal)[signal]
  allBS.bg <- sum(sets[which(!(sets$feature %in% significant)),"sites"])
  allBS.sig <- sum(sets[which(sets$feature %in% significant),"sites"])
  res <- t(vapply(split(sets,sets$set), set1=significant, bs.sig=allBS.sig, 
                  bs.bg=allBS.bg, FUN.VALUE=numeric(8), 
                  FUN=function(x,set1,bs.sig,bs.bg){
                    x <- as.data.frame(x)
                    w <- which(as.character(x$feature) %in% set1)
                    xin <- sum(x[w,"sites"])
                    xout <- sum(x[which(!(as.character(x$feature) %in% set1)),"sites"])
                    mm <- matrix(round(c(xin,bs.sig-xin,xout,bs.bg-xout)),nrow=2)
                    p1 <- fisher.test(mm,alternative="greater")$p.value[[1]]
                    p2 <- fisher.test(mm,alternative="less")$p.value[[1]]
                    c(  overlap=length(unique(x[w,"feature"])),
                        BS.in=xin,
                        otherBS.in=bs.sig-xin,
                        BS.inBG=xout,
                        enrichment=log2(((1+xin)/(bs.sig-xin))/((1+xout)/(bs.bg-xout))),
                        otherBS.inBG=bs.bg-xout,
                        under.pvalue=p2,
                        pvalue=p1
                    )
                  }))
  res[res[,"pvalue"]==0,"pvalue"] <- min(res[res[,"pvalue"]>0,"pvalue"])/10
  res <- as.data.frame(res)
  res <- res[order(res$pvalue),]
  for(i in 1:4) res[[i]] <- as.integer(res[[i]])
  res$FDR <- p.adjust(res$pvalue,method="fdr")
  res
}






#' siteoverlap
#'
#' Applies Fisher's test to the number of binding sites among a set of features 
#' (vs all other binding sites in that set of features).
#' Of note, this application violates the assumptions of Fisher's exact test.
#'
#' @param signal A named logical vector indicating which features are 
#' differentially-expressed
#' @param sets A data.frame with at least the following columns: 
#' 'set', 'feature', 'sites'.
#' @param alternative Test alternative (defaults to 'greater' to test for
#' over-representation)
#'
#' @return a data.frame.
#'
#' @export
siteoverlap <- function(signal, sets,
                        alternative=c("greater","less","two.sided")){
  alternative <- match.arg(alternative)
  if(!.checkSets(sets, "sites", matrixAlternative="numeric")){
    sets <- .aggregateByFamilies(sets)
    sets <- .setsToScoreMatrix(signal, sets, column="sites", keepSparse=TRUE)
  }else{
    sets <- sets[row.names(sets) %in% names(signal),,drop=FALSE]
  }
  signal <- signal[row.names(sets)]
  BS.in <- colSums(sets[which(signal),,drop=FALSE])
  otherBS.in <- sum(sets[which(signal),,drop=FALSE])-BS.in
  BS.inBG <- colSums(sets[which(!signal),,drop=FALSE])
  otherBS.inBG <- sum(sets[which(!signal),,drop=FALSE])-BS.inBG
  enrichment <- ((1+BS.in)/(otherBS.in))/((1+BS.inBG)/(otherBS.inBG))
  expected <- sum(signal)*(BS.in+BS.inBG)/nrow(sets)
  res <- data.frame( overlap=colSums(sets[which(signal),,drop=FALSE]>0),
                     sites.overlap=BS.in,
                     #BS.inBG=BS.inBG,
                     #otherBS.in=otherBS.in,
                     enrichment=round(log2(enrichment),3),
                     pvalue=fisher.test.p(BS.in, otherBS.in,
                                          BS.inBG, otherBS.inBG,
                                          alternative=alternative) )
  res[res$pvalue==0,"pvalue"] <- min(res[res$pvalue>0,"pvalue"])/10
  res$FDR <- p.adjust(res$pvalue, method="fdr")
  res[order(res$FDR,res$pvalue),]
}




#' plMod
#'
#' Fits linear models between a signal and set variables. This is a flexible but
#' slow implementations, we recommend \link{ebayes} for a faster one.
#'
#' @param signature A named numeric vector
#' @param sets A data.frame with at least the following columns: 'set', 
#' 'feature', and the variable specified by `var`.
#' @param var The independent variable to use, either 'sites' or 'score'.
#' @param correctForLength Logical; whether to correct for length (i.e. total
#'  sites). Defaults to TRUE if `var='sites'`, FALSE otherwise.
#'
#' @return a data.frame
#' @importFrom stats .lm.fit p.adjust
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
  sets$sets <- as.character(sets$set)
  res <- vapply(split(as.data.frame(sets),sets$set), fcs=signature, 
                  FUN.VALUE=numeric(2), FUN=function(x,fcs, minSize){
    x <- x[!duplicated(x$feature),]
    row.names(x) <- x$feature
    x2 <- x[names(fcs),var]
    x2[which(is.na(x2))] <- 0
    if(!is.null(cfl)){
      x2 <- cbind(x2=x2, cfl=cfl)
    }else{
      x2 <- as.matrix(x2)
    }
    mod <- try(.lm.fit(x2, fcs), silent=TRUE)
    if(is(mod,"try-error")) return(rep(NA_real_,2))
    c(mod$coefficients[1], tryCatch(.lm.pval(mod)[1], error=function(e) NA))
  })
  res <- as.data.frame(t(res))
  colnames(res) <- c("coefficient","pvalue")
  res$FDR <- p.adjust(res$pvalue)
  res[order(res$FDR,res$pvalue),]
}

#' @export
#' @rdname plMod
modsites <- function(x, sets, correctForLength=TRUE, ...){
  plMod(x, sets, correctForLength=correctForLength, ...)
}

#' @export
#' @rdname plMod
modscore <- function(x, sets, ...) plMod(x, sets, var="score", ...)


#' woverlap
#'
#' Weighted hypergeometric test adjusted for total number of miRNA bindings 
#' sites. This is done through the `goseq` package.
#'
#' @param signal A named logical vector indicating the features of interest 
#' (e.g. which features are differentially-expressed)
#' @param sets A data.frame with at least the following columns: 
#' 'set', 'feature'.
#' @param method Method for handling then length bias, default "Wallenius".
#' See `?goseq` for more detail.
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
    ag <- rowsum(sets$sites, sets$feature)
    ag <- ag[intersect(row.names(ag),names(bd)),]
    bd[names(ag)] <- ag
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

#' ks
#'
#' enrichment analysis using a Kolmogorov-Smirnov test on the signal.
#'
#' @param signal A named logical vector indicating the features of interest 
#' (e.g. which features are differentially-expressed)
#' @param sets A data.frame with at least the following columns: 'set', 'feature'.
#'
#' @return a data.frame.
#' @export
#' @importFrom stats ks.test
ks <- function(signal, sets){
  res <- vapply(split(sets$feature,sets$set), FUN.VALUE=numeric(1), FUN=function(x){
    x <- intersect(x, names(signal))
    if(length(x)==0 || length(x)==length(signal)) return(NA_real_)
    ks <- suppressWarnings(try(ks.test(signal[x], signal[setdiff(names(signal),x)])$p.value, silent=T))
    if(is(ks,"try-error")) return(NA_real_)
    ks
  })
  res <- data.frame(row.names = names(res), pvalue=as.numeric(res))
  res$FDR <- p.adjust(res$pvalue,method="fdr")
  res[order(res$FDR,res$pvalue),]
}

#' mw
#'
#' miRNA targets enrichment analysis using a Mann-Whitney / Wilcoxon test on 
#' the targets' foldchanges.
#'
#' @param signal A named logical vector indicating the features of interest 
#' (e.g. which features are differentially-expressed)
#' @param sets A data.frame with at least the following columns: 
#' 'set', 'feature'.
#'
#' @return a data.frame.
#' @export
#' @importFrom stats wilcox.test
mw <- function(signal, sets){
  res <- vapply(split(sets$feature,sets$set), FUN.VALUE=numeric(1), FUN=function(x){
    x <- intersect(x, names(signal))
    ks <- suppressWarnings(try(wilcox.test(signal[x], signal[setdiff(names(signal),x)])$p.value, silent=T))
    if(is(ks,"try-error")) return(NA)
    ks
  })
  res <- data.frame(row.names = names(res), pvalue=as.numeric(res))
  res$FDR <- p.adjust(res$pvalue,method="fdr")
  res[order(res$FDR,res$pvalue),]
}


#' gsea
#'
#' @param signal A named numeric vector indicating the features' response (e.g. 
#' logFC)
#' @param sets A data.frame with at least the following columns:
#' 'set', 'feature'.
#' @param maxSize The maximum number of elements in a set to be tested (default 
#' 500). If the test takes too long to run, consider setting this.
#' @param nperm The number of permutations, default 2000. The more permutations,
#'  the more precise the p-values, in particular for small p-values.
#'
#' @return a data.frame.
#'
#' @export
# @importFrom fgsea fgseaMultilevel
gsea <- function(signal, sets, maxSize=3000, nperm=2000, ...){
  library(fgsea)
  sets <- lapply(split(sets$feature,sets$set), tested=names(signal), 
                 FUN=function(x,tested){ intersect(unique(x),tested) })
  res <- fgseaMultilevel(sets, signal, nPermSimple=nperm, minSize=4, 
                         maxSize=maxSize)
  res <- as.data.frame(res[order(res$padj,res$pval),])
  colnames(res)[1:5] <- c("family","pvalue","FDR","ES","normalizedEnrichment")
  colnames(res)[8] <- "features"
  row.names(res) <- res[,1]
  return(res[,-1])
}


.censorScore <- function(x){
  x <- 0.1-x
  if(length(w <- which(x>0.9))>0)
    x[w] <- 0.9+ecdf(x[w])(x[w])/10
  x
}

.TS2regulon <- function(x, likelihood="score"){
  x <- as.data.frame(x)
  if(likelihood %in% colnames(x)){
    if(likelihood=="score"){
      x$likelihood <- .censorScore(x[[likelihood]])
    }else{
      x$likelihood <- x[[likelihood]]
    }
  }else{
    x$likelihood <- 1L
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
#' @importFrom viper msviper
#' @export
areamir <- function(signal, sets, ...){
  vi <- viper::msviper(signal, regulon=.TS2regulon(as.data.frame(sets)), ...,
                       verbose=FALSE)
  vi2 <- as.data.frame(vi$es[c("nes","p.value")])
  colnames(vi2) <- c("enrichment","pvalue")
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
#' 'set', 'feature', and (if binary=FALSE) 'score'. Alternatively, a sparse
#' logical (if binary=TRUE) or numeric matrix with features as rows and sets
#' as columns.
#' @param binary Logical; whether to consider target prediction as binary.
#' @param alpha elastic net mixing param (0=ridge, 1=lasso)
#' @param do.plot Logical, whether to plot coefficients against lambda
#' @param use.intercept Logical, whether to use an intercept in the model.
#' @param keepAll Logical, whether to return all families.
#'
#' @return A DataFrame.
#' @import sparseMatrixStats
#' @export
# @importFrom zetadiv glm.cons
regmir <- function(signal, sets, binary=NULL, alpha=1, do.plot=FALSE, 
                   use.intercept=FALSE, keepAll=TRUE){
  if(is.null(names(signal))) stop("`signal` should be a named vector!")
  suppressPackageStartupMessages({
    library(glmnet)
    library(zetadiv)
  })
  if(is.null(binary)){
    if(is.data.frame(sets) || is(sets,"DFrame")){
      binary <- !("score" %in% colnames(sets))
    }else{
      binary <- is.logical(sets) || is(sets,"lgCMatrix")
    }
  }
  
  # prepare the target matrix
  if(binary){
    if(!.checkSets(sets, matrixAlternative="logical")){
      sets <- .setsToScoreMatrix(signal, sets, column=NULL, keepSparse=TRUE)
    }
  }else{
    sets <- .setsToScoreMatrix(signal, sets, keepSparse=TRUE)
  }
  sets <- sets[row.names(sets) %in% names(signal),]
  sets <- .reduceBm(sets)
  signal <- signal[row.names(sets)]
  
  # regularized regression with cross-validation
  if(isLogistic <- is.logical(signal)){
    fits <- cv.glmnet(sets, signal, standardize=FALSE, alpha=alpha, 
                      family="binomial", lower.limits=0, intercept=use.intercept)
  }else{
    fits <- cv.glmnet(sets, signal, standardize=FALSE, alpha=alpha, 
                      family="gaussian", intercept=use.intercept)
  }
  
  if(do.plot){
    layout(matrix(1:2,nrow=1))
    plot(fits)
    plot(fits$glmnet.fit, label=TRUE)
  }
  
  # we extract the miRNAs selected by the best most regularized glmnet fit:
  co <- .decideLambda(fits)

  # p <- fixedLassoInf(sets, signal, co[,1], lambda=lambda, intercept=use.intercept, 
  #                    family=ifelse(isLogistic,"binomial","gaussian"), alpha=alpha-0.01)
  # res <- data.frame(row.names=names(p$vars), beta=co[names(p$vars),1], pvalue=p$pv)
  # print(res)
  
  beta <- co[,1]
  names(beta) <- row.names(co)
  signald <- data.frame( y=signal, as.matrix(sets[,names(beta),drop=FALSE]) )
  if(!binary){
    signald$median <- sparseMatrixStats::rowMedians(sets)
    if(all(signald$median==0)) signald$median <- rowMeans(sets)
  } 
  
  # new fit to get significance estimates
  if(use.intercept){
    form <- y~.
  }else{
    form <- y~0+.
  }
  res <- tryCatch({
    if(isLogistic){
      mod <- glm.cons( form, data=signald, family="binomial", cons=1, cons.inter=-1)
    }else{
      mod <- lm( form, data=signald )
    }
    
    # we extract the coefficients and p-values, and reorganize the output:
    res <- coef(summary(mod))
    res[order(res[,4]),,drop=FALSE]
  }, error=function(e){
    warning(e)
    data.frame(row.names=names(beta), beta=beta, z=rep(NA_real_, length(beta)), pvalue=1)
  })
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
  
  res <- res[row.names(res)!="median",]
  # res <- DataFrame(res)
  # we adjust using all features as number of comparisons
  if(nrow(res)>0){
    res$FDR <- p.adjust(res$pvalue, n=max(nrow(res),ncol(sets),na.rm=TRUE))
    # res$features <- CharacterList(lapply(split(TS$feature, TS$family)[row.names(res)], 
    #                                 y=names(signal)[signal], FUN=intersect))
  }
  res
}


#' @export
#' @rdname regmir
regmir.bb <- function(signal, sets, ...) regmir(signal, sets, binary=TRUE, ...)
#' @export
#' @rdname regmir
regmir.bc <- function(signal, sets, ...) regmir(signal, sets, binary=FALSE, ...)
#' @export
#' @rdname regmir
regmir.cc <- function(signal, sets, ...) regmir(signal, sets, binary=FALSE, ...)

.decideLambda <- function(fits){
  l <- fits$lambda.1se
  if(fits$nzero[fits$lambda==l]<2){
    w1 <- which(fits$lambda==fits$lambda.min)
    w2 <- which(fits$lambda==fits$lambda.1se)
    l <- fits$lambda[floor((w1-w2)/2+w2)]
  }
  if(fits$nzero[fits$lambda==l]<2) l <- fits$lambda.min
  co <- coef(fits, s=l)
  co[co[,1]!=0,,drop=FALSE]
}

.lm.pval <- function(m){
  rdf <- length(m$residuals) - ncol(m$qr)
  x <- sum(m$residuals^2)/rdf
  se <- sqrt(diag(chol2inv(m$qr))*x)
  2*pt(abs(m$coef/se),rdf,lower.tail=FALSE)
}

# checks for duplicated bm columns for regmir
.reduceBm <- function(bm){
  # split columns into groups having the same colSums
  si <- split(seq_len(ncol(bm)), colSums(bm))
  if(!any(lengths(si)>1)) return(bm)
  ll <- lapply(si[lengths(si)>1], FUN=function(i){
    # for each group of >1 columns:
    # calculate pairwise distances
    d <- as.matrix(dist(t(bm[,i])))
    # remove redundant sets
    d <- d[!duplicated(d),,drop=FALSE]
    # extract sets
    lapply(seq_len(nrow(d)), FUN=function(x){
      colnames(d)[which(d[x,]==0)]
    })
  })
  ll <- unlist(ll, recursive=FALSE)
  # remove things that aren't duplicated
  ll <- ll[lengths(ll)>1]
  if(length(ll)==0) return(bm)
  old.names <- sapply(ll, FUN=function(x) x[1]) # those columns we'll rename
  toRemove <- sapply(ll, FUN=function(x) x[-1]) # those columns we remove
  bm <- bm[,setdiff(colnames(bm), unlist(toRemove))]
  bm2 <- bm[,old.names]
  colnames(bm2) <- sapply(ll, FUN=function(x) paste(sort(x),collapse="."))
  cbind(bm[,setdiff(colnames(bm), old.names)],bm2)
}

#' @import Matrix
.setsToScoreMatrix <- function(signal, sets, column="score", keepSparse=FALSE){
  if(!is.factor(sets$feature)) sets$feature <- as.factor(sets$feature)
  signal <- signal[names(signal) %in% levels(sets$feature)]
  sets$set <- droplevels(as.factor(sets$set))
  if(is.null(column)){
    column <- TRUE
  }else{
    column <- sets[[column]]
  }
  bm <- sparseMatrix(i=as.integer(sets$feature), j=as.integer(sets$set), 
                     x=column, dim=c(length(levels(sets$feature)), 
                                     length(levels(sets$set))),
                     dimnames=list(levels(sets$feature),levels(sets$set)))
  bm <- bm[names(signal),,drop=FALSE]
  if(!keepSparse) bm <- as.matrix(bm)
  bm
}

#' ebayes
#'
#' A wrapper around \code{\link[limma]{lmFit}} and \code{\link[limma]{eBayes}}
#' for fitting the feature scores of many sets against a signal.
#' 
#' @param signal A named numeric vector indicating the signal (e.g. logFC) of
#' each feature.
#' @param sets A data.frame of annotation, with at least the following columns: 
#' `set`, `feature` and `score`. Alternatively, a sparse numeric matrix, with sets as 
#' columns and features as rows (dimensions must be named).
#' @param use.intercept Logical; whether to use an intercept in the model.
#' 
#' @return A data.frame.
#'
#' @importFrom limma lmFit topTable eBayes
#' @importFrom stats model.matrix
#' @import sparseMatrixStats
#' @export
ebayes <- function(signal, sets, use.intercept=FALSE){
  if(!.checkSets(sets, "score", matrixAlternative="numeric"))
    sets <- .setsToScoreMatrix(signal, sets, keepSparse=TRUE)
  signal <- signal[names(signal) %in% row.names(sets)]
  sets <- sets[names(signal),]
  meds <- sparseMatrixStats::rowMedians(sets)
  if(all(meds==0)) meds <- rowMeans(sets)
  if(use.intercept){
    mm <- model.matrix(~meds+signal)
  }else{
    mm <- model.matrix(~0+meds+signal)
  }
  fit <- lmFit(as.matrix(t(sets)), mm)
  res <- as.data.frame(topTable(eBayes(fit), coef="signal", number = Inf)[,c(1,4,5)])
  colnames(res) <- c("coefficient", "pvalue", "FDR")
  res
}

#' lmadd
#'
#' Complementary linear fits between signal and sets' feature scores
#' 
#' @param signal A named numeric vector indicating the signal (e.g. logFC) of
#' each feature.
#' @param sets A data.frame of annotation, with at least the following columns: 
#' `set`, `feature` and `score`. Alternatively, a sparse numeric matrix, with sets as 
#' columns and features as rows (dimensions must be named).
#' @param use.intercept Logical; whether to use an intercept in the model.
#' @param calc.threshold Minimum p-value threshold for the individual fit in 
#' order to calculate a complementary coefficient.
#' @param comb.threshold Minimum p-value threshold for the individual fit in 
#' order to include as covariate for other coefficients.
#' 
#' @return A data.frame.
#'
#' @importFrom limma lmFit topTable eBayes
#' @import sparseMatrixStats
#' @importFrom stats coef lm.fit
#' @export
lmadd <- function(signal, sets, use.intercept=FALSE, calc.threshold=0.2, 
                  comb.threshold=0.05){
  if(!.checkSets(sets, "score", matrixAlternative="numeric"))
    sets <- .setsToScoreMatrix(signal, sets, keepSparse=TRUE)
  signal <- signal[names(signal) %in% row.names(sets)]
  sets <- sets[names(signal),]
  meds <- sparseMatrixStats::rowMedians(sets)
  if( skipMeds <- (sum(meds==0)/length(meds))>0.9){
    if(use.intercept){
      mm <- model.matrix(~signal)
    }else{
      mm <- model.matrix(~0+signal)
    }
  }else{
    if(use.intercept){
      mm <- model.matrix(~meds+signal)
    }else{
      mm <- model.matrix(~0+meds+signal)
    }
  }
  fit <- eBayes(lmFit(as.matrix(t(sets)), mm))
  res1 <- as.data.frame(topTable(fit, coef=1+!skipMeds+use.intercept, 
                                 number=Inf)[,c(1,1,4,4)])
  colnames(res1) <- c("independent.coef", "combined.coef", 
                      "independent.pvalue", "combined.pvalue")
  res1$combined.coef <- res1$independent.coef <- 0
  res1$combined.pvalue[-1] <- 1
  i <- 1
  p <- res1$independent.pvalue[1]
  sets <- data.frame(signal=signal, intercept=1, median=meds, as.matrix(sets), 
                     check.names=FALSE)
  if(skipMeds) sets$median <- NULL
  while(i<=nrow(res1) && p<=calc.threshold){
    feats <- row.names(res1)[i]
    if(use.intercept) feats <- c("intercept", feats)
    fit <- lm.fit(as.matrix(sets[,feats,drop=FALSE]), signal)
    res1$independent.coef[i] <- rev(fit$coefficients)[1]
    feats <- c("signal",row.names(res1)[seq_len(i)])
    if(i==1 && !skipMeds) feats <- c("median",feats)
    if(use.intercept) feats <- c("intercept",feats)
    fit <- lm(signal~.,data=sets[,feats])
    co <- coef(summary(fit))
    res1$combined.pvalue[i] <- p <- rev(co[,4])[1]
    if(p<comb.threshold){
      res1$combined.coef[rev(seq_len(i))] <- rev(co[,1])[seq_len(i)]
    }else{
      res1$combined.coef[i] <- rev(co[,1])[1]
    }
    i <- i+1
  }
  res1$FDR <- p.adjust(res1$combined.pvalue)
  res1
}

#' fisher.test.p
#' 
#' Fast p-values from multiple Fisher's exact tests
#'
#' @param a,b,c,f Vectors containing the four positions of the 2x2 matrix
#' @param alternative greater, less, or 'two.sided' (default)
#'
#' @return A vector of p-values
#' @importFrom stats dhyper
#' @export
fisher.test.p <- function (a, b, c, d, 
                           alternative=c("two.sided", "greater", "less")){
  fn <- switch( match.arg(alternative), 
     less = function(x,m,n,k) phyper(x, m, n, k), 
     greater = function(x,m,n,k) phyper(x - 1, m, n, k, lower.tail=FALSE), 
     two.sided = function(x,m,n,k){
       lo <- max(0, k - n)
       support <- seq(from=lo, to=min(k, m))
       d <- dhyper(support, m, n, k, log = TRUE)
       d <- exp(d - max(d))
       d <- d/sum(d)
       sum(d[d <= d[x - lo + 1] * (1 + 10^(-7))])
     }
  )
  mapply(FUN=fn, a, a+c, b+d, a+b)
}

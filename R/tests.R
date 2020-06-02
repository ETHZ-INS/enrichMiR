#' EA
#'
#' Traditional overlap (Fisher test) between differentially-expressed features and miRNA targets.
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#' @param testOnlyAnnotated Whether to excluded features that are bound by no miRNA (default FALSE).
#'
#' @return a data.frame.
#'
#' @export
EA <- function(signal,TS, minSize=5, testOnlyAnnotated=FALSE){
  library(data.table)
  tested <- names(signal)
  significant <- names(signal)[signal]
  if(testOnlyAnnotated) tested <- intersect(tested,unique(as.character(TS$feature)))
  significant <- intersect(significant,tested)
  sets <- split(TS,TS$family)
  sets <- sets[lengths(sets)>=minSize]
  res <- lapply(sets, set1=significant, universe=tested, FUN=function(x,set1,universe){
    set2 <- unique(as.character(x$feature))
    set2 <- intersect(set2,universe)
    ov <- intersect(set2,set1)
    expected <- length(set1)*length(set2)/length(universe)
    list(  family=as.character(x$family[1]),
           annotated=length(set2),
           overlap=length(ov),
           expected=round(expected,2),
           enrichment=round(log2((length(ov)+0.1)/expected),2),
           under.pvalue=.overlap.prob(set1,set2,universe,lower=T),
           over.pvalue=.overlap.prob(set1,set2,universe),
           features=ov
    )
  })
  feats <- CharacterList(lapply(res, FUN=function(x) x$features))
  res <- DataFrame(rbindlist(lapply(res, FUN=function(x) x[-length(x)])))
  res$FDR <- p.adjust(res$over.pvalue,method="fdr")
  res$features <- feats
  row.names(res) <- res$family
  res[order(res$FDR,res$over.pvalue),-1]
}


#' plMod
#'
#' @param dea A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#' @param var The independent variable to use, either 'sites' or 'score'.
#' @param correctForLength Logical; whether to correct for UTR size (using total bindings as a proxy). 
#' Defaults to TRUE if `var='sites'`, FALSE otherwise.
#'
#' @return a data.frame
#' @export
plMod <- function(dea, TS, minSize=5, var="sites", correctForLength=(var=="sites")){
  library(MASS)
  fcs <- dea$logFC
  names(fcs) <- row.names(dea)
  TS <- aggregate(TS,by=list(family=TS$family, feature=TS$feature),FUN=function(x){ if(is.numeric(x)) return(max(x,na.rm=T)); x[[1]] })
  if(correctForLength){
    ag <- aggregate(TS$sites,by=list(feature=TS$feature),FUN=sum)
    cfl <- ag[,2]
    names(cfl) <- ag[,1]
    cfl <- cfl[names(fcs)]
    cfl[which(is.na(cfl))] <- 0
    names(cfl) <- names(fcs)
  }else{
    cfl <- NULL
  }
  res <- t(sapply(split(TS,TS$family,drop=F), fcs=fcs, minSize=minSize, cfl=cfl, FUN=function(x,fcs, minSize, cfl){
    r <- c(x$family[1],x$rep.miRNA[1],NA,NA)
    if(nrow(x)<minSize) return(r)
    x <- x[!duplicated(x),]
    row.names(x) <- x$feature
    x2 <- x[names(fcs),var]
    x2[which(is.na(x2))] <- 0
    if(is.null(cfl)){
      mod <- try(rlm(fcs~x2+0),silent=T)
    }else{
      mod <- try(rlm(fcs~x2+cfl+0),silent=T)
    }
    if(!is(mod,"try-error")) r[3:4] <- c(coef(mod)["x2"], summary(aov(mod))[[1]]["x2","Pr(>F)"])
    return(r)
  }))
  res <- DataFrame(res)
  colnames(res) <- c("family","rep.miRNA","logFC","pvalue")
  res$logFC <- as.numeric(as.character(res$logFC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$FDR <- p.adjust(res$pvalue)
  row.names(res) <- res$family
  res[order(res$FDR,res$pvalue),-1]
}

#' wEA
#'
#' Weighted hypergeometric test adjusted for total number of miRNA bindings sites. (This is done through the `goseq` package)
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#' @param testOnlyAnnotated Whether to excluded features that are bound by no miRNA (default FALSE).
#' @param method Method for handling then length bias, default "Wallenius". See `?goseq` for more detail.
#'
#' @return a data.frame.
#'
#' @importFrom goseq nullp goseq
#' @export
wEA <- function(signal,TS, minSize=5, testOnlyAnnotated=FALSE, method="Wallenius"){
  ag <- aggregate(TS$sites,by=list(gene=TS$feature),FUN=sum)
  row.names(ag) <- ag$gene
  if(testOnlyAnnotated) signal <- signal[intersect(names(signal), ag$gene)]
  ag <- ag[names(signal),]
  bd <- ag$x
  names(bd) <- ag$gene
  np <- nullp(signal, bias.data = bd, plot.fit=FALSE)
  g2c <- TS[,c("feature","family")]

  res <- goseq(np, gene2cat=g2c, method=method, use_genes_without_cat=!testOnlyAnnotated)
  colnames(res) <- c("family","over.pvalue","under.pvalue","overlap","numInCat")

  # TEMPORARY enrichment value
  significant <- names(signal)[signal]
  res$enrichment <- round(log2(res$overlap/(length(significant)*(res$numInCat/length(signal)))),2)
  
  res$FDR <- p.adjust(res$over.pvalue, method="fdr")
  colnames(res)[5] <- "annotated"
  ll <- split(TS$feature,TS$family)
  res <- DataFrame(res)
  res$features <- CharacterList(lapply( ll[as.character(res$family)],
                                        y=significant, FUN=function(x,y){
                                                            sort(intersect(x,y))
                                                          }))
  row.names(res) <- res$family
  res[,-1]
}

#' michael
#'
#' Applies Fisher's test to the number of miRNA binding sites among a set of features (vs all other binding sites in that set of features).
#' Of note, this application violates the assumptions of Fisher's exact test.
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param significant The set of differentially-expressed features.
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#' @param testOnlyAnnotated Whether to excluded features that are bound by no miRNA (default FALSE).
#' @param alternative Alternative tested, default "greater" overlap.
#'
#' @return a data.frame.
#'
#' @export
michael <- function(signal, TS, minSize=3, testOnlyAnnotated=FALSE){
  library(data.table)
  tested <- names(signal)
  significant <- names(signal)[signal]
  if(testOnlyAnnotated) tested <- intersect(tested,unique(as.character(TS$feature)))
  TS <- TS[which(TS$feature %in% tested),]
  significant <- intersect(significant,tested)
  allBS.bg <- sum(TS[which(!(TS$feature %in% significant)),"sites"])
  allBS.sig <- sum(TS[which(TS$feature %in% significant),"sites"])
  ll <- split(TS,TS$family,drop=F)
  res <- lapply(ll, set1=significant, bs.sig=allBS.sig, bs.bg=allBS.bg, alternative=alternative, FUN=function(x,set1,bs.sig,bs.bg,alternative){
    w <- which(as.character(x$feature) %in% set1)
    xin <- sum(x[w,"sites"])
    xout <- sum(x[which(!(as.character(x$feature) %in% set1)),"sites"])
    mm <- matrix(c(xin,bs.sig-xin,xout,bs.bg-xout),nrow=2)
    p1 <- fisher.test(mm,alternative="greater")$p.value[[1]]
    p2 <- fisher.test(mm,alternative="less")$p.value[[1]]
    list(   family=as.character(x$family[1]),
            annotated=nrow(x),
            overlap=length(unique(x[w,"feature"])),
            BS.in=xin,
            otherBS.in=bs.sig-xin,
            BS.inBG=xout,
            enrichment=log2((xin/(bs.sig-xin))/(xout/(bs.bg-xout))),
            otherBS.inBG=bs.bg-xout,
            under.pvalue=p2,
            over.pvalue=p1
    )
  })
  res <- res[which(sapply(res,length)>0)]
  res <- DataFrame(rbindlist(res))
  res <- res[which(res[,"annotated"]>=minSize),]
  res$fisher.pvalue <- as.numeric(unlist(res$over.pvalue))
  res$FDR <- p.adjust(res$over.pvalue,method="fdr")
  res$features <- CharacterList(lapply( ll[as.character(res$family)],
                                        y=significant, FUN=function(x,y){
                                          sort(intersect(x,y))
                                        }))
  row.names(res) <- res$family
  res[order(res$FDR,res$over.pvalue),-1]
}

#' KS
#'
#' miRNA targets enrichment analysis using a Kolmogorov-Smirnov test on the targets' foldchanges. As all alternatives are considered, significance might not always be consistent with differential miRNA activity.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#'
#' @return a data.frame.
#'
#' @export
KS <- function(DEA, TS, minSize=5){
  library(data.table)
  res <- lapply(split(TS,TS$family), DEA=DEA, FUN=function(x,DEA){
    set <- unique(as.character(x$feature))
    set <- intersect(set,row.names(DEA))
    ks <- suppressWarnings(try(ks.test(DEA[set,"logFC"], DEA[setdiff(row.names(DEA),set),"logFC"])$p.value, silent=T))
    if(is(ks,"try-error")) ks <- NA
    list(family=as.character(x$family[1]),
         annotated=length(set),
         ks.pvalue=ks)
  })
  res <- DataFrame(rbindlist(res))
  res <- res[which(res[,"annotated"]>=minSize),]
  res$FDR <- p.adjust(res$ks.pvalue,method="fdr")
  row.names(res) <- res$family
  res[order(res$FDR),-1]
}

#' KS2
#'
#' miRNA targets enrichment analysis using a Kolmogorov-Smirnov test on the targets' foldchanges, treating upregulated and downregulated genes separately, and applying a one-sided test to each. Significance here is more likely to incidate an effect consistent with miRNA activity than in a normal KS-test.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#'
#' @return a data.frame.
#'
#' @export
KS2 <- function(DEA, TS, minSize=5){
  library(data.table)
  res <- lapply(split(TS,TS$family), DEA=DEA, FUN=function(x,DEA){
    set <- unique(as.character(x$feature))
    d1 <- DEA[which(DEA$logFC>0),]
    d2 <- DEA[which(DEA$logFC<0),]
    set1 <- intersect(set,row.names(d1))
    set2 <- intersect(set,row.names(d2))
    ks1 <- suppressWarnings(try(ks.test(d1[set1,"logFC"], DEA[setdiff(row.names(d1),set1),"logFC"],alternative="less")$p.value, silent=T))
    ks2 <- suppressWarnings(try(ks.test(d2[set2,"logFC"], DEA[setdiff(row.names(d2),set2),"logFC"],alternative="greater")$p.value, silent=T))
    if(is(ks1,"try-error")) ks1 <- NA
    if(is(ks2,"try-error")) ks2 <- NA
    list(family=as.character(x$family[1]),
         annotated=length(set),
         ks.pvalue.down=ks2,
         ks.pvalue.up=ks1
    )
  })
  res <- DataFrame(rbindlist(res))
  if(!(nrow(res)>1)) return(res)
  res <- res[which(res[,"annotated"]>=minSize),]
  res$ks.pvalue.down <- unlist(res$ks.pvalue.down)
  res$ks.pvalue.up <- unlist(res$ks.pvalue.up)
  res$FDR <- apply(matrix(p.adjust(as.numeric(as.matrix(res[,c("ks.pvalue.down","ks.pvalue.up")])),method="fdr"),ncol=2),1,FUN=min)
  row.names(res) <- res$family
  res[order(res$FDR,apply(res[,grep("pvalue",colnames(res))],1,FUN=min)),-1]
}


#' MW
#'
#' miRNA targets enrichment analysis using a Mann-Whitney / Wilcoxon test on the targets' foldchanges.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#'
#' @return a data.frame.
#'
#' @export
MW <- function(DEA, TS, minSize=5){
  library(data.table)
  res <-lapply(split(TS,TS$family), DEA=DEA, FUN=function(x,DEA){
    set <- unique(as.character(x$feature))
    set <- intersect(set,row.names(DEA))
    mw <- try(wilcox.test(DEA[set,"logFC"], DEA[setdiff(row.names(DEA),set),"logFC"])$p.value, silent=T)
    if(is(mw,"try-error")) mw <- NA
    list(family=as.character(x$family[1]),
         annotated=length(set),
         wilcox.pvalue=mw)
  })
  res <- DataFrame(rbindlist(res))
  res <- res[which(res[,"annotated"]>=minSize),]
  res$FDR <- p.adjust(res$wilcox.pvalue,method="fdr")
  row.names(res) <- res$family
  res[order(res$FDR),-1]
}


#' gsea
#'
#' miRNA target enrichment analysis among differentially-expressed genes using GSEA.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of elements in a set to be tested (default 5).
#' @param maxSize The maximum number of elements in a set to be tested (default 500). If the test takes too long to run, consider setting this.
#' @param fdr.thres The FDR threshold below which genes are considered; default 0.2.
#' @param nperm The number of permutations, default 2000. The more permutations, the more precise the p-values, in particular for small p-values.
#'
#' @return a data.frame.
#'
#' @importFrom fgsea fgsea
#' @export
gsea <- function(DEA, TS, minSize=5, maxSize=300, fdr.thres=0.5, nperm=2000){
  DEA <- DEA[which(DEA$FDR<=fdr.thres & !is.na(DEA$logFC)),]
  w <- which(is.infinite(DEA$logFC))
  if(length(w)>0) DEA$logFC[w] <- max(abs(DEA$logFC[-w]))*sign(DEA$logFC[w])
  fcs <- DEA$logFC
  names(fcs) <- row.names(DEA)
  sets <- lapply(split(TS$feature,TS$family),tested=names(fcs),FUN=function(x,tested){ intersect(unique(x),tested) })
  sets <- sets[which(sapply(sets,length)>=minSize)]
  res <- fgsea(sets, fcs, nperm, minSize=minSize, maxSize=maxSize)
  res <- DataFrame(res[order(res$padj,res$pval),])
  colnames(res)[1:5] <- c("family","pvalue","FDR","ES","normalizedEnrichment")
  colnames(res)[8] <- "features"
  row.names(res) <- res$family
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
  
  res <- sapply(split(TS[,c("family","sites")],TS$feature,drop=F), fams=unique(TS$family), bgs=bgs, FUN=function(x, fams, bgs){
    p <- rep(1,length(fams))
    names(p) <- fams
    tbs.in <- sum(x$sites)
    tbs <- sum(bgs)
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
  w <- which(x>0.9)
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
  lapply(split(x,x$family, drop=TRUE), FUN=function(x){
    y <- list(  tfmode=rep(-1,nrow(x)),
                likelihood=x$likelihood )
    lapply(y, FUN=function(a){
      names(a) <- x$feature
      a
    })	
  })
}

#' aREAmir
#'
#' analytic Rank-based Enrichment Analysis using a conversion of targetScan 
#' scores as weights.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of elements in a set to be tested (default 5).
#' @param maxSize The maximum number of elements in a set to be tested (default 500). If the test takes too long to run, consider setting this.
#' @param fdr.thres The FDR threshold below which genes are considered; default 0.2.
#' @param nperm The number of permutations, default 2000. The more permutations, the more precise the p-values, in particular for small p-values.
#'
#' @return a data.frame.
#'
#' @export
aREAmir <- function(dea, TS, minSize=5, pleiotropy=FALSE){
  sig <- -log10(dea$FDR)*sign(dea$logFC)
  names(sig) <- row.names(dea)
  vi <- viper::msviper(sig, regulon=TS2regulon(as.data.frame(TS)), minsize=minSize, pleiotropy=pleiotropy, verbose=FALSE)
  vi2 <- DataFrame(vi$es[c("nes","size","p.value","nes.bt")])
  colnames(vi2)[3] <- "pvalue"
  vi2$FDR <- p.adjust(vi2$pvalue)
  vi2 <- vi2[order(vi2$pvalue),]
  #vi2$miRNAs <- sapply(row.names(vi2), fam=metadata(TS)$families, FUN=function(x, fam) names(fam)[which(fam==x)])
  vi2
}

#' regmir
#'
#' miRNA enrichment analysis using regularized regression
#'
#' @param signal A vector of logical or numeric values, with gene symbols as names.
#' @param TS The targetscan target predictions object
#' @param binary Logical; whether to consider target prediction as binary (default). The
#' alternative is to use the score.
#' @param alpha elastic net mixing param (0=ridge, 1=lasso)
#' @param do.plot Logical, whether to plot coefficients against lambda
#' @param use.intercept Logical, whether to use an intercept in the model.
#'
#' @return A DataFrame.
#' @importFrom IRanges CharacterList
#' @import glmnet zetadiv S4Vectors 
#' @export
regmir <- function(signal, TS, binary=TRUE, alpha=1, do.plot=FALSE, use.intercept=FALSE){
  if(is.null(names(signal))) stop("`signal` should be a named vector!")
  
  # prepare the target matrix
  if(binary){
    bm <- sapply(split(TS$feature, TS$family), FUN=function(x) names(signal) %in% x)
  }else{
    TS2 <- aggregate(TS[,"score",drop=FALSE], by=as.data.frame(TS[,c("family","feature")]), FUN=min)
    TS2 <- split(TS2[,c("score","feature")], TS2$family)
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
  res <- res[grep("^\\(Intercept\\)$|FALSE$", row.names(res), invert=TRUE),,drop=FALSE]
  row.names(res) <- gsub("TRUE","",row.names(res))
  res <- DataFrame(res[order(res[,4]),,drop=FALSE])
  colnames(res) <- c("beta","stderr",ifelse(isLogistic,"z","t"),"pvalue")
  # we adjust using all features as number of comparisons
  res$FDR <- p.adjust(res$pvalue, n=ncol(bm))
  res$features <- CharacterList(lapply(split(TS$feature, TS$family)[row.names(res)], 
                                    y=names(signal)[signal], FUN=intersect))
  res
}

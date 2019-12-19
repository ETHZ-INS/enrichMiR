#' EA
#'
#' Traditional overlap (Fisher test) between differentially-expressed features and miRNA targets.
#'
#' @param tested The set of tested features (e.g. genes)
#' @param significant The set of differentially-expressed features.
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#' @param testOnlyAnnotated Whether to excluded features that are bound by no miRNA (default FALSE).
#'
#' @return a data.frame.
#'
#' @export
EA <- function(tested,significant,TS, minSize=5, testOnlyAnnotated=FALSE){
  library(data.table)
  tested <- as.character(tested)
  significant <- as.character(significant)
  if(testOnlyAnnotated) tested <- intersect(tested,unique(as.character(TS$feature)))
  significant <- intersect(significant,tested)
  res <- lapply(split(TS,TS$family), set1=significant, universe=tested, FUN=function(x,set1,universe){
    set2 <- unique(as.character(x$feature))
    set2 <- intersect(set2,universe)
    ov <- intersect(set2,set1)
    expected <- length(set1)*length(set2)/length(universe)
    list(  family=as.character(x$family[1]),
           annotated=length(set2),
           overlap=length(ov),
           expected=round(expected,2),
           enrichment=round(length(ov)/expected,2),
           under.pvalue=.overlap.prob(set1,set2,universe,lower=T),
           over.pvalue=.overlap.prob(set1,set2,universe),
           features=paste(ov,collapse=", ")
    )
  })
  res <- as.data.frame(rbindlist(res))
  res <- res[which(res[,"annotated"]>=minSize),]
  res$FDR <- p.adjust(res$over.pvalue,method="fdr")
  res <- res[,c(1:(ncol(res)-2),ncol(res),ncol(res)-1)]
  res[order(res$FDR,res$over.pvalue),]
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
  res <- as.data.frame(res)
  colnames(res) <- c("family","rep.miRNA","logFC","pvalue")
  res$logFC <- as.numeric(as.character(res$logFC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$FDR <- p.adjust(res$pvalue)
  res[order(res$FDR,res$pvalue),]
}

#' wEA
#'
#' Weighted hypergeometric test adjusted for total number of miRNA bindings sites. (This is done through the `goseq` package)
#'
#' @param tested The set of tested features (e.g. genes)
#' @param significant The set of differentially-expressed features.
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#' @param testOnlyAnnotated Whether to excluded features that are bound by no miRNA (default FALSE).
#' @param method Method for handling then length bias, default "Wallenius". See `?goseq` for more detail.
#'
#' @return a data.frame.
#'
#' @export
wEA <- function(tested,significant,TS, minSize=5, testOnlyAnnotated=FALSE, method="Wallenius"){
  library(goseq)
  ag <- aggregate(TS$sites,by=list(gene=TS$feature),FUN=sum)
  row.names(ag) <- ag$gene
  if(testOnlyAnnotated) tested <- intersect(tested, ag$gene)
  ag <- ag[tested,]
  bd <- ag$x
  names(bd) <- ag$gene
  de <- tested %in% significant
  names(de) <- tested
  np <- nullp(de, bias.data = bd, plot.fit=FALSE)
  g2c <- TS[,c("feature","family")]

  res <- goseq(np, gene2cat=g2c, method=method, use_genes_without_cat=!testOnlyAnnotated)
  colnames(res) <- c("family","over.pvalue","under.pvalue","overlap","numInCat")

  # TEMPORARY enrichment value
  res$enrichment <- round(res$overlap/(length(significant)*(res$numInCat/length(tested))),2)
  
  res$FDR <- p.adjust(res$over.pvalue, method="fdr")
  colnames(res)[5] <- "annotated"
  ll <- split(TS$feature,TS$family)
  res$features <- sapply(ll[as.character(res$family)],y=significant,FUN=function(x,y){ paste(sort(intersect(x,y)),collapse=", ") })
  res
}

#' michael
#'
#' Applies Fisher's test to the number of miRNA binding sites among a set of features (vs all other binding sites in that set of features).
#' Of note, this application violates the assumptions of Fisher's exact test.
#'
#' @param tested The set of tested features (e.g. genes)
#' @param significant The set of differentially-expressed features.
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#' @param testOnlyAnnotated Whether to excluded features that are bound by no miRNA (default FALSE).
#' @param alternative Alternative tested, default "greater" overlap.
#'
#' @return a data.frame.
#'
#' @export
michael <- function(tested, significant, TS, minSize=3, testOnlyAnnotated=FALSE){
  library(data.table)
  tested <- as.character(tested)
  significant <- as.character(significant)
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
            enrichment=(xin/(bs.sig-xin))/(xout/(bs.bg-xout)),
            otherBS.inBG=bs.bg-xout,
            under.pvalue=p2,
            over.pvalue=p1
    )
  })
  res <- res[which(sapply(res,length)>0)]
  res <- as.data.frame(rbindlist(res))
  res <- res[which(res[,"annotated"]>=minSize),]
  res$fisher.pvalue <- as.numeric(unlist(res$over.pvalue))
  res$FDR <- p.adjust(res$over.pvalue,method="fdr")
    ll <- split(TS$feature,TS$family)
  res$features <- sapply(ll[as.character(res$family)],y=significant,FUN=function(x,y){ paste(sort(intersect(x,y)),collapse=", ") })
  
  res[order(res$FDR,res$over.pvalue),]
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
  res <- as.data.frame(rbindlist(res))
  res <- res[which(res[,"annotated"]>=minSize),]
  res$FDR <- p.adjust(res$ks.pvalue,method="fdr")
  res[order(res$FDR),]
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
  res <- as.data.frame(rbindlist(res))
  if(!(nrow(res)>1)) return(res)
  res <- res[which(res[,"annotated"]>=minSize),]
  res$ks.pvalue.down <- unlist(res$ks.pvalue.down)
  res$ks.pvalue.up <- unlist(res$ks.pvalue.up)
  res$FDR <- apply(matrix(p.adjust(as.numeric(as.matrix(res[,c("ks.pvalue.down","ks.pvalue.up")])),method="fdr"),ncol=2),1,FUN=min)
  res[order(res$FDR,apply(res[,grep("pvalue",colnames(res))],1,FUN=min)),]
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
  res <- as.data.frame(rbindlist(res))
  res <- res[which(res[,"annotated"]>=minSize),]
  res$FDR <- p.adjust(res$wilcox.pvalue,method="fdr")
  res[order(res$FDR),]
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
#' @export
gsea <- function(DEA, TS, minSize=5, maxSize=300, fdr.thres=0.5, nperm=2000){
  library(fgsea)
  DEA <- DEA[which(DEA$FDR<=fdr.thres & !is.na(DEA$logFC)),]
  w <- which(is.infinite(DEA$logFC))
  if(length(w)>0) DEA$logFC[w] <- max(abs(DEA$logFC[-w]))*sign(DEA$logFC[w])
  fcs <- DEA$logFC
  names(fcs) <- row.names(DEA)
  sets <- lapply(split(TS$feature,TS$family),tested=names(fcs),FUN=function(x,tested){ intersect(unique(x),tested) })
  sets <- sets[which(sapply(sets,length)>=minSize)]
  res <- fgsea(sets, fcs, nperm, minSize=minSize, maxSize=maxSize)
  res <- res[order(res$padj,res$pval),]
  colnames(res)[1:5] <- c("family","pvalue","FDR","ES","normalizedEnrichment")
  res$leadingEdge <- sapply(res$leadingEdge, FUN=function(x){ paste(x,collapse=", ") })
  return(as.data.frame(res))
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
#' @export
geneBasedTest <- function(features, TS){
  library(aggregation)
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

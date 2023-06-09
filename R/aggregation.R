#' Aggregate the various tests
#'
#' Significance is assessed with a permutation of the p percentiles.
#'
#' @param er An `enrich.results` object
#' @param niter The number of iterations for the non-parametric p-value estimation
#' @param discordantFactor The factor by which to increase the mean p quantile of
#'   discordant enrichments
#'
#' @return A results dataframe
#' @export
aggregateTests <- function(er, niter=20000, discordantFactor=5){
  stopifnot(is(er,"enrich.results"))
  tests <- names(er@res)
  res <- getResults(er)
  tcon <- tests[grep("up|down",tests,inver=TRUE)]
  tbin <- tests[grep("up|down",tests)]
  tbin <- unique(gsub("\\.up|\\.down","",tbin))
  if(!(length(tcon)>0 & length(tbin)>0))
    stop("This function requries at least one binary and one continuous tests,",
         " ideally more.")
  colnames(res) <- gsub("coefficient","enrichment",colnames(res))
  colnames(res) <- gsub("\\.over","",colnames(res))
  res$minFDR <- apply(res[,grep("FDR",colnames(res))],1,FUN=min)
  for(f in grep("pvalue$",colnames(res),value=TRUE)){
    x <- res[[f]]
    res[[gsub("pvalue","pq",f)]] <- ecdf(x)(x)
  }

  wM <- function(p, e, pq=FALSE){
    if(pq){
      w <- apply(1-as.matrix(p),2,FUN=function(x) ecdf(x)(x))
    }else{
      w <- 1-as.matrix(p)
    }
    w <- w/rowSums(w,na.rm=TRUE)
    rowSums(apply(as.matrix(e),2,FUN=function(x){
      x/sd(x,na.rm=TRUE)
    })*w, na.rm=TRUE)
  }

  res$enr.c <- wM(res[,paste0(tcon,".pq"),drop=FALSE],
                  res[,paste0(tcon,".enrichment"),drop=FALSE])

  for(f in tbin){
    res[[paste0(f,".pq")]] <- ifelse(res$enr.c>0,res[[paste0(f,".up.pq")]],
                                     res[[paste0(f,".down.pq")]])
    res[[paste0(f,".enrichment")]] <- ifelse(res$enr.c>0,res[[paste0(f,".up.enrichment")]],
                                             res[[paste0(f,".down.enrichment")]])
  }
  res <- res[,grep("\\.up|\\.down|\\.over|\\.under",colnames(res),invert=TRUE)]
  res$enr.b <- wM(res[,paste0(tbin,".pq"),drop=FALSE],
                  res[,paste0(tbin,".enrichment"),drop=FALSE])
  res$enr <- rowMeans(cbind(res$enr.b,res$enr.c))

  pq <- grep("pq$",colnames(res),value=TRUE)
  res$meanPQ <- rowMeans(as.matrix(res[,pq]),na.rm=TRUE)
  # penalize inconsistent enrichments
  w <- which(sign(res$enr.b)!=sign(res$enr.c))
  res$meanPQ[w] <- res$meanPQ[w]*discordantFactor
  xnull <- as.data.frame(lapply(res[,pq], FUN=function(j){
    x <- sample(j, length(j)*niter, replace=TRUE)
  }))
  xnullfn <- ecdf(rowMeans(as.matrix(xnull),na.rm=TRUE))
  res$pvalue <- xnullfn(res$meanPQ)
  res$FDR <- p.adjust(res$pvalue)
  res$enrichment <- res$enr
  cols <- c("members","set_size","enrichment","pvalue","FDR.geomean","FDR.mean","FDR")
  res <- res[,intersect(cols,colnames(res))]
  res[order(res$pvalue),]
}

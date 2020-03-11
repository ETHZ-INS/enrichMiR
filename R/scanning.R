#' findSeedMatches
#'
#' @param seqs A character vector of sequences in which to look. If DNA, will be
#' complemented before matching.
#' @param seeds A character vector of 7-nt seeds to look for. If RNA, will be 
#' reversed and complemented before matching. If DNA, they are assumed to be
#' the target sequence to look for. Alternatively, a list of objects of class
#' `KdModel` can be given.
#' @param shadow Integer giving the shadow, i.e. the number of nucleotides
#'  hidden at the beginning of the sequence (default 0)
#' @param keepMatchSeq Logical; whether to keep the sequence (including flanking
#' dinucleotides) for each seed match (default FALSE).
#' @param seedtype Either RNA, DNA or 'auto' (default)
#' @param BP Pass `BiocParallel::MulticoreParam(ncores)` to enable 
#' multithreading.
#'
#' @return A GRanges of all matches
#' 
#' @import BiocParallel Biostrings GenomicRanges
#' @importFrom stringr str_locate_all str_sub
#' @export
#'
#' @examples
#' # we create mock RNA sequences and seeds:
#' seqs <- sapply(1:10, FUN=function(x) paste(sample(strsplit("ACGU", "")[[1]], 
#'                                      1000, replace=TRUE),collapse=""))
#' names(seqs) <- paste0("seq",1:length(seqs))
#' seeds <- c("AAACCAC", "AAACCUU")
#' m <- findSeedMatches(seqs, seeds)
  findSeedMatches <- function( seqs, seeds, shadow=0, keepMatchSeq=FALSE,
                             seedtype=c("auto", "RNA","DNA"), BP=NULL){
  library(GenomicRanges)
  library(stringr)
  library(BiocParallel)
  library(Biostrings)
  if(is.null(names(seqs))) names(seqs) <- paste0("seq",seq_along(seqs)) 
  seedtype <- match.arg(seedtype)
  seqtype <- .guessSeqType(seqs)
  if(is.list(seeds) && all(sapply(seeds, FUN=function(x) is(x,"KdModel")))){
    if(is.null(names(seeds)))
      stop("If `seeds` is a list of kd models, it should be named.")
    if(seedtype=="RNA" || seqtype=="RNA") 
      stop("If `seeds` is a list of kd models, both the seeds and the target
sequences should be in DNA format.")
  }else{
    if(is.null(names(seeds))) names(seeds) <- seeds
    if(seedtype=="auto") seedtype <- .guessSeqType(seeds)
  }
  n <- names(seeds)
  if(seedtype=="RNA"){
    message("Matching reverse complements of the seeds...")
    seeds <- as.character(reverseComplement(RNAStringSet(seeds)))
  }else{
    message("Matching the given seeds directly...")
  }
  if(seqtype=="RNA"){
    seeds <- gsub("T", "U", seeds)
  }else{
    seeds <- gsub("U", "T", seeds)
  }
  names(seeds) <- n
  if(is.null(BP)) BP <- SerialParam()
  if(shadow>0) seqs <- substr(seqs, shadow+1, sapply(seqs, nchr))
  seqs <- seqs[sapply(seqs,nchar)>=min(sapply(seeds,nchar))]
  m <- lapply(seeds, seqs=seqs, FUN=function(seed,seqs){
    mod <- NULL
    if(is(seed,"KdModel")){
      mod <- seed
      seed <- substring(seed$xlevels$sr,2)
    }
    seed <- paste0(".?.?.?",substr(seed,2,7),".?.?.?")
    pos <- stringr::str_locate_all(seqs, seed)
    if(sum(sapply(pos,nrow))==0) return(GRanges())
    y <- GRanges( rep(names(seqs), sapply(pos,nrow)), 
                  IRanges( start=unlist(lapply(pos,FUN=function(x) x[,1])),
                           end=unlist(lapply(pos,FUN=function(x) x[,2])) ) )
    y$sequence <- factor(unlist(mapply(x=seqs, pos=pos, FUN=function(x,pos){
      stringr::str_sub(x, pos[,1], pos[,2])
    })))
    y <- characterizeSeedMatches( y,
                                  ifelse(is.null(mod),seed,mod$canonical.seed),
                                  mod )
    if(!keepMatchSeq) y$sequence <- NULL
    y
  })
  for(s in names(seeds)) m[[s]]$seed <- s
  m <- sort(unlist(GRangesList(m)))
  m$seed <- factor(m$seed)
  m$type <- factor(m$type)
  if(shadow>0){
    end(m) <- end(m)+shadow
    start(m) <- start(m)+shadow
  }
  m
}

.guessSeqType <- function(x, use.subset=TRUE){
  seqs <- x[sample.int(length(x),min(length(x),10))]
  u <- any(grepl("U",seqs))
  t <- any(grepl("T",seqs))
  if(t && u) stop("Sequences contain both T and U!")
  if(t || u) return(ifelse(u,"RNA","DNA"))
  if(length(seqs)>10) return(.guessSeqType(x, FALSE))
  warning("Sequences contain neither T or U - assuming they are DNA...")
  return("DNA")
}

#' characterizeSeedMatches
#'
#' @param x A factor or character vector of matches, or a data.frame or GRanges
#' containing a `sequence` column.
#' @param seed The seed which sequences of `x` should match.
#' @param mod An optional object of class 'KdModel'.
#'
#' @return A data.frame, or a `GRanges` if `x` is a `GRanges`
#' @export
#'
#' @examples
#' characterizeSeedMatches(c("UAAACCACCC","CGAACCACUG"), "AAACCAC")
characterizeSeedMatches <- function(x, seed=NULL, kd.model=NULL){
  if(is.null(seed)){
    if(is.null(kd.model))
      stop("At least one of `seed` or `kd.model` must be given.")
    seed <- kd.model$canonical.seed
  }
  seed <- as.character(seed)
  if(length(seed)>1 || nchar(seed)!=7) 
    stop("`seed` must contain a single string of 7 characters.")
  if(!is.character(x) && !is.factor(x)){
    y <- characterizeSeedMatches(x$sequence, seed, kd.model)
    if(is(x,"GRanges")){
      x@elementMetadata <- cbind(x@elementMetadata, y)
    }else{
      x <- cbind(x,y)
    }
    return(x)
  }
  if(any(duplicated(x))){
    # we resolve type for unique sequences
    if(is.factor(x)){
      y <- characterizeSeedMatches(levels(x), seed, kd.model)
    }else{
      y <- characterizeSeedMatches(unique(x), seed, kd.model)
    }
    y <- y[as.character(x),,drop=FALSE]
    row.names(y) <- NULL
    return(y)
  }
  d <- data.frame( row.names=x, 
                   type=sapply(as.character(x), seed=seed, FUN=.getMatchType) )
  if(!is.null(kd.model)) d$log_kd <- predictKD(row.names(d), kd.model)
  d
}

.getMatchType <- function(x, seed){
  if(grepl(paste0(seed,"A"),x,fixed=TRUE)) return("8mer")
  if(grepl(seed,x,fixed=TRUE)) return("7mer-m8")
  seed6 <- substr(seed,2,7)
  if(grepl(paste0(seed6,"A"),x,fixed=TRUE)) return("7mer-A1")
  if(grepl(seed6,x,fixed=TRUE)) return("6mer")
  "non-canonical"
}


#' getKdModel
#' 
#' Summarizes the binding affinity of 12-mers using linear models.
#'
#' @param kd A data.frame with at least the columns "X12mer" and "log_kd"
#'
#' @return A linear model of class `KdModel`
#' @export
getKdModel <- function(kd){
  if(is.character(kd) && length(kd)==1) kd <- read.delim(kd, header=TRUE)
  if("mirseq" %in% colnames(kd)){
    mirseq <- as.character(kd$mirseq[1])
    #seed <- substr(mirseq, 2,8)
    seed <- as.character(complement(DNAString(substr(mirseq, 2,8)))
    kd <- kd[,c("X12mer","log_kd")]
  }else{
    mirseq <- seed <- NULL
  }
  kd <- kd[grep("X",kd$X12mer,invert=TRUE),]
  pwm <- Biostrings::consensusMatrix(
    as.character(rep(kd$X12mer, floor( (10^(-kd$log_kd))/10 ))),
    as.prob=TRUE
  )
  fields <- c("sr","A","fl")
  if(!all(fields %in% colnames(kd)))
    kd <- prep12mers(kd)
  fields <- c(fields, "log_kd")
  if(!all(fields %in% colnames(kd))) stop("Malformed `kd` data.frame.")
  mod <- lm( log_kd~sr*A+fl, data=kd, model=FALSE, weights=(-kd$log_kd)^2, 
             x=FALSE, y=FALSE )
  mod$residuals <- NULL
  mod$fitted.values <- NULL
  mod$weights <- NULL
  mod$assign <- NULL
  mod$effects <- NULL
  mod$qr <- list(pivot=mod$qr$pivot)
  mod$mirseq <- mirseq
  mod$canonical.seed <- seed
  mod$pwm <- pwm
  class(mod) <- c("KdModel", class(mod))
  mod
}

#' prep12mers
#'
#' @param x A vector of 12-mers, or a data.frame containing at least the columns
#' "X12mer" and "log_kd"
#' @param mod An optional linear model summarizing the kd activity.
#' @param maxSeedMedian Max median log_kd for alternative seed inclusion.
#'
#' @return A data.frame
#' @export
prep12mers <- function(x, mod=NULL, maxSeedMedian=-1.2){
  if(is.data.frame(x)){
    if(!all(c("X12mer","log_kd") %in% colnames(x)))
      stop("`x` should be a character vector or a data.frame with the columns ",
           "'X12mer' and 'log_kd'")
    x <- x[grep("X",x$X12mer,invert=TRUE),]
    x <- cbind(x[,"log_kd",drop=FALSE], prep12mers(x$X12mer))
    ag <- aggregate(x$log_kd, by=list(seed=x$sr), FUN=median)
    ag <- as.character(ag[ag$x <= maxSeedMedian,1])
    return( rbind( data.frame( log_kd=0, sr="other", A=FALSE, 
                               fl=levels(x$fl)[1] ),
                   x[x$sr %in% ag,] ) )
  }
  x <- as.character(x)
  sr <- sapply(x, FUN=function(x) substr(x, 3,9))
  if(!is.null(mod)){
    if(is.null(mod$xlevels$sr)) stop("The model contains no seed levels.")
    sr.lvls <- mod$xlevels$sr
    sr[!(sr %in% sr.lvls)] <- "other"
  }else{
    sr.lvls <- c("other",unique(sr))
  }
  sr <- factor(sr, sr.lvls)
  fl <- sapply(x, FUN=function(x) paste0(substr(x, 1, 2),substr(x, 11, 12)))
  data.frame(sr=sr, A=sapply(x, FUN=function(x) substr(x, 10, 10)=="A"),
             fl=factor(fl, levels=getKmers(4)), row.names=NULL)
}

#' getKmers
#'
#' Returns all combinations of `n` elements of `from`
#'
#' @param n Number of elements
#' @param from Letters sampled
#'
#' @return A character vector
#' @export
#'
#' @examples
#' getKmers(3)
getKmers <- function(n=4, from=c("G", "C", "T", "A")){
  apply(expand.grid(lapply(seq_len(n), FUN=function(x) from)),
        1,collapse="",FUN=paste)
}

#' predictKD
#'
#' @param kmer The 12-mer sequences for which affinity should be predicted
#' @param mod  A `KdModel`
#'
#' @return A vector of the same length as `kmer` indicating the corresponding
#' predicted affinities.
#' @export
predictKD <- function(kmer, mod){
  if(!is(mod,"KdModel")) stop("`mod` should be of class `KdModel`.")
  kd <- suppressWarnings(predict(mod, prep12mers(kmer, mod)))
  kd[which(kd>0)] <- 0
  kd
}

#' plotKdModel
#'
#' @param mod A `KdModel`
#' @param what Either 'seeds', 'logo', or 'both' (default)
#'
#' @return A plot
#' @export
plotKdModel <- function(mod, what=c("both","seeds","logo")){
  library(ggplot2)
  what <- match.arg(what)
  if(what=="seeds"){
    co <- -coefficients(mod)[paste0("sr",mod$xlevels$sr[-1])]
    names(co) <- gsub("^sr","",names(co))
    co <- sort(co)
    co <- data.frame(seed=factor(names(co), names(co)), log_kd=as.numeric(co))
    co$type <- sapply(as.character(co$seed), seed=mod$canonical.seed, .getMatchType)
    return( ggplot(co, aes(seed, log_kd, fill=type)) + geom_col() + 
              coord_flip() + ylab("-log_kd") )
  }
  if(what=="logo") return(seqLogo::seqLogo(mod$pwm))
  library(cowplot)
  plot_grid( plotKdModel(mod, "seeds"),
             grid::grid.grabExpr(plotKdModel(mod, "logo")),
             nrow=2)
}
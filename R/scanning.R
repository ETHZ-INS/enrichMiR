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
  }
  if(is.null(BP)) BP <- SerialParam()
  if(shadow>0) seqs <- substr(seqs, shadow+1, sapply(seqs, nchar))
  seqs <- seqs[sapply(seqs,nchar)>=min(sapply(seeds,nchar))]
  seqnms <- factor(names(seqs), names(seqs))
  seqs <- paste0("xxx",seqs,"xxx")
  names(seqs) <- levels(seqnms)
  m <- bplapply(seeds, seqs=seqs, BPPARAM=BP, FUN=function(seed,seqs){
    if(is(seed,"KdModel")){
      seed2 <- substring(seed$xlevels$sr[-1],2)
    }else{
      seed2 <- substr(seed,2,7)
    }
    # look-around matching to get overlapping seeds
    pos <- gregexpr( paste0("(?=",paste(unique(seed2),collapse="|"),")"),
                     seqs, perl=TRUE )
    pos <- lapply(lapply(pos, as.numeric), y=-1, setdiff)
    if(sum(sapply(pos,length))==0) return(NULL)
    GRanges( rep(seqnms, sapply(pos,length)), 
             IRanges( start=unlist(pos), width=6 ) )
  })
  m <- m[!sapply(m,is.null)]
  mseed <- factor(rep(names(m),sapply(m,length)))
  m <- unlist(GRangesList(m))
  m$seed <- mseed
  m <- unlist(GRangesList(lapply(split(m,seqnames(m)), FUN=function(r){
    r$sequence <- stringr::str_sub( seqs[[as.numeric(seqnames(r[1]))]], 
                                    start(r)-3, end(r)+3 )
    r
  })))
  row.names(m) <- NULL
  start(m) <- start(m)-3
  end(m) <- end(m)+3
  m <- unlist(GRangesList(bplapply(split(m, m$seed), BPPARAM=BP, FUN=function(x){
    seed <- seeds[[as.character(x$seed[1])]]
    mod <- NULL
    if(is(seed,"KdModel")){
      mod <- seed
      seed <- seed$canonical.seed
    }
    x <- characterizeSeedMatches( x, seed, mod)
    # for overlapping seeds, keep only the best one
    x <- x[order(x$log_kd),]
    .removeOverlapping(x)
  })))
  if(!keepMatchSeq) m$sequence <- NULL
  m <- sort(m)
  m$type <- factor(m$type)
  if(shadow>0){
    end(m) <- end(m)+shadow
    start(m) <- start(m)+shadow
  }
  m
}

# keeps the first
.removeOverlapping <- function(x){
  red <- reduce(x, with.revmap=TRUE)
  red <- red[lengths(red$revmap)>1]
  if(length(red)==0) return(x)
  il <- lengths(red$revmap)
  red <- unlist(red$revmap)
  toKeep <- red[cumsum(il)-il+1]
  x[-setdiff(red, toKeep)]
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
  if(grepl(paste0("[ACGT]",seed6,"A"),x)) return("7mer-A1")
  if(grepl(paste0("[ACGT]",seed6),x)) return("6mer")
  if(grepl(seed6,x,fixed=TRUE)) return("offset 6mer")
  "non-canonical"
}


#' getKdModel
#' 
#' Summarizes the binding affinity of 12-mers using linear models.
#'
#' @param kd A data.frame with at least the columns "X12mer" and "log_kd", or the path to
#' such a data.frame
#' @param name The optional name of the miRNA
#'
#' @return A linear model of class `KdModel`
#' @export
getKdModel <- function(kd, name=NULL){
  if(is.character(kd) && length(kd)==1){
    if(is.null(name)) name <- gsub("\\.rds$","",basename(kd),ignore.case=TRUE)
    kd <- read.delim(kd, header=TRUE)
  }
  if("mirseq" %in% colnames(kd)){
    mirseq <- as.character(kd$mirseq[1])
    seed <- as.character(reverseComplement(DNAString(substr(mirseq, 2,8))))
    name <- as.character(kd$mir[1])
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
  w <- (1-kd$log_kd)^2
  w[which(w<0.5)] <- 0.5
  mod <- lm( log_kd~sr*A+fl, data=kd, model=FALSE, weights=w, x=FALSE, y=FALSE )
  mod$cor.with.cnn <- cor(mod$fitted.values, kd$log_kd)
  mod$mae.with.cnn <- median(abs(mod$fitted.values-kd$log_kd))
  mod$residuals <- NULL
  mod$fitted.values <- NULL
  mod$weights <- NULL
  mod$assign <- NULL
  mod$effects <- NULL
  mod$qr <- list(pivot=mod$qr$pivot)
  mod$name <- name
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
  d <- data.frame( sr=sr, A=sapply(x, FUN=function(x) substr(x, 10, 10)=="A"),
                   fl=factor(fl, levels=getKmers(4)), row.names=NULL )
  d$A[is.na(d$A)] <- FALSE
  fl <- sort(coef(mod)[grep("fl",names(coef(mod)))])
  d$fl[is.na(d$fl)] <- names(fl)[floor(length(fl)/2)]
  d
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
    coe <- coefficients(mod)
    medfl <- median(coe[grep("^fl",names(coe),value=TRUE)],na.rm=TRUE)
    co <- -coe["(Intercept)"]-coe[paste0("sr",mod$xlevels$sr[-1])]-medfl
    names(co) <- gsub("^sr","",names(co))
    co <- sort(co)
    co <- data.frame(seed=factor(names(co), names(co)), log_kd=as.numeric(co))
    co$type <- sapply(as.character(co$seed), seed=mod$canonical.seed, .getMatchType)
    coA <- co
    aint <- coe[paste0("sr",coA$seed,":ATRUE")]
    aint[is.na(aint)] <- 0
    coA$log_kd <- - coe["ATRUE"] - aint
    coA$type <- "+A"
    co <- rbind(co,coA)
    co$type <- factor(co$type, c("+A","7mer-m8","6mer","offset 6mer","non-canonical"))
    # mer8 <- co[nrow(co),,drop=FALSE]
    # mer8$log_kd <- mer8$log_kd-coe[paste0("sr",mer8$seed,":ATRUE")]-coe["ATRUE"]
    # mer8$type="8mer"
    # co <- rbind(co, mer8)
    # co$type <- factor(co$type, c("8mer","7mer-m8","6mer","non-canonical"))
    p <- ggplot(co, aes(seed, log_kd, fill=type)) + geom_col() + 
      coord_flip() + ylab("-log_kd")
    if(!is.null(mod$name)) p <- p + ggtitle(mod$name)
    return( p )
  }
  if(what=="logo") return(seqLogo::seqLogo(mod$pwm))
  library(cowplot)
  plot_grid( plotKdModel(mod, "seeds"),
             grid::grid.grabExpr(plotKdModel(mod, "logo")),
             nrow=2)
}

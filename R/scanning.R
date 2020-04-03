#' findSeedMatches
#'
#' @param seqs A character vector of sequences in which to look. If DNA, will be
#' complemented before matching.
#' @param seeds A character vector of 7-nt seeds to look for. If RNA, will be 
#' reversed and complemented before matching. If DNA, they are assumed to be
#' the target sequence to look for. Alternatively, a list of objects of class
#' `KdModel` or an object of class `CompressedKdModelList` can be given.
#' @param shadow Integer giving the shadow, i.e. the number of nucleotides
#'  hidden at the beginning of the sequence (default 0)
#' @param keepMatchSeq Logical; whether to keep the sequence (including flanking
#' dinucleotides) for each seed match (default FALSE).
#' @param minDist Integer specifying the minimum distance between matches of the same 
#' miRNA (default 1). Closer matches will be reduced to the highest-affinity of the two.
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
findSeedMatches <- function( seqs, seeds, shadow=0, keepMatchSeq=FALSE, minDist=1,
                             seedtype=c("auto", "RNA","DNA"), BP=NULL){
  library(GenomicRanges)
  library(stringr)
  library(BiocParallel)
  library(Biostrings)
  if(is.null(names(seqs))) names(seqs) <- paste0("seq",seq_along(seqs)) 
  seedtype <- match.arg(seedtype)
  seqtype <- .guessSeqType(seqs)
  if(is(seeds, "CompressedKdModelList")) seeds <- decompressKdModList(seeds)
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
  m <- unlist(GRangesList(bplapply( split(m,seqnames(m),drop=TRUE), 
                                    BPPARAM=BP, FUN=function(r){
    if(length(r)>0)
      r$sequence <- stringr::str_sub( seqs[[as.numeric(seqnames(r[1]))]], 
                                      start(r)-3, end(r)+3 )
    r
  })))
  row.names(m) <- NULL
  start(m) <- start(m)-3
  end(m) <- end(m)-3
  m <- unlist(GRangesList(bplapply(split(m, m$seed), BPPARAM=BP, FUN=function(x){
    seed <- seeds[[as.character(x$seed[1])]]
    if(is(seed,"KdModel")){
      x <- characterizeSeedMatches( x, seed$canonical.seed, seed)
      x <- x[order(x$log_kd),]
    }else{
      x <- characterizeSeedMatches( x, seed)
      x <- x[order(x$type),]
    }
    # for overlapping seeds, keep only the best one
    .removeOverlapping(x, minDist=minDist)
  })))
  if(!keepMatchSeq) m$sequence <- NULL
  m <- sort(m)
  m$type <- factor(m$type)
  if(shadow>0){
    end(m) <- end(m)+shadow
    start(m) <- start(m)+shadow
  }
  names(m) <- NULL
  m$log_kd <- round(m$log_kd, 3)
  m
}

# keeps the first
.removeOverlapping <- function(x, minDist=1L){
  red <- reduce(x, with.revmap=TRUE, min.gapwidth=minDist)
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
  d$type <- factor(d$type, 
                   c("8mer","7mer-m8","7mer-A1","6mer","offset 6mer","non-canonical"))
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

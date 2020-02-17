#' findSeedMatches
#'
#' @param seqs A character vector of sequences in which to look. If DNA, will be
#' complemented before matching.
#' @param seeds A character vector of 7-nt seeds to look for. If RNA, will be 
#' reversed and complemented before matching. If DNA, they are assumed to be
#' the target sequence to look for.
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
  if(is.null(names(seeds))){ n <- seeds }else{ n <- names(seeds) }
  seqtype <- .guessSeqType(seqs)
  if(seedtype=="auto") seedtype <- .guessSeqType(seeds)
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
    pos <- stringr::str_locate_all(seqs, paste0(".?.?",substr(seed,2,7),".?.?"))
    if(sum(sapply(pos,nrow))==0) return(GRanges())
    y <- GRanges( rep(names(seqs), sapply(pos,nrow)), 
                  IRanges( start=unlist(lapply(pos,FUN=function(x) x[,1])),
                           end=unlist(lapply(pos,FUN=function(x) x[,2])) ) )
    y$sequence <- factor(unlist(mapply(x=seqs, pos=pos, FUN=function(x,pos){
      stringr::str_sub(x, pos[,1], pos[,2])
    })))
    y <- characterizeSeedMatches(y, seed)
    if(!keepMatchSeq) y$sequence <- NULL
    y
  })
  for(s in names(seeds)) m[[s]]$seed <- seeds[[s]]
  m <- sort(unlist(GRangesList(m)))
  m$seed <- factor(m$seed)
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
#' @param seed The seed which sequences of `x` match.
#' @param iScore The base scores attributed to different types of matches,
#' before considering dinucleotides. Should be a numeric vector with the
#' following names: "8mer", "7mer-m8", "7mer-A1", "6mer"
#'
#' @return A data.frame, or a `GRanges` if `x` is a `GRanges`
#' @export
#'
#' @examples
#' characterizeSeedMatches(c("UAAACCACCC","CGAACCACUG"), "AAACCAC")
characterizeSeedMatches <- function(x, seed, iScore=c("8mer"=2.5, "7mer-m8"=1.8,
                                                      "7mer-A1"=1.5, "6mer"=1)){
  seed <- as.character(seed)
  if(length(seed)>1 || nchar(seed)!=7) 
    stop("`seed` must contain a single string of 7 characters.")
  if(!is.character(x) && !is.factor(x)){
    y <- characterizeSeedMatches(x$sequence, seed)
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
      y <- characterizeSeedMatches(levels(x), seed)
    }else{
      y <- characterizeSeedMatches(unique(x), seed)
    }
    y <- y[as.character(x),,drop=FALSE]
    row.names(y) <- NULL
    return(y)
  }
  seed.length <- length(seed)
  d <- as.data.frame(t(sapply(as.character(x), FUN=function(x){
    seed6 <- substr(seed,2,7)
    if(grepl(paste0(seed,"A"),x,fixed=TRUE)){
      type <- "8mer"
    }else if(grepl(seed,x,fixed=TRUE)){
      type <- "7mer-m8"
    }else if(grepl(paste0(seed6,"A"),x,fixed=TRUE)){
      type <- "7mer-A1"
    }else{
      type <- "6mer"
    }
    score <- as.numeric(iScore[type])
    flanking <- strsplit(x, paste0(substr(seed,1,1),"?",seed6))[[1]]
    c(type=type, score=score)
  })))
  d$score <- as.numeric(as.character(d$score))
  row.names(d) <- x
  d
}

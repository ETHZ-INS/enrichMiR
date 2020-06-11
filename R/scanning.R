#' findSeedMatches
#'
#' @param seqs A character vector or `XStringSet`` of sequences in which to look.
#' @param seeds A character vector of 7-nt seeds to look for. If RNA, will be 
#' reversed and complemented before matching. If DNA, they are assumed to be
#' the target sequence to look for. Alternatively, a list of objects of class
#' `KdModel` or an object of class `CompressedKdModelList` can be given.
#' @param seedtype Either RNA, DNA or 'auto' (default)
#' @param shadow Integer giving the shadow, i.e. the number of nucleotides
#'  hidden at the beginning of the sequence (default 0)
#' @param minLogKd Minimum log_kd value to keep (default 1). Set to 0 to disable.
#' @param types The type of sites to return (default all)
#' @param max.noncanonical.motifs The maximum number of non-canonical motifs to search 
#' for (default all)
#' @param keepMatchSeq Logical; whether to keep the sequence (including flanking
#' dinucleotides) for each seed match (default FALSE).
#' @param minDist Integer specifying the minimum distance between matches of the same 
#' miRNA (default 1). Closer matches will be reduced to the highest-affinity. To 
#' disable the removal of overlapping features, use `minDist=-Inf`.
#' @param fastRemoveOverlaps Whether to use the fast (rather than exact) method to remove 
#' overlaps. Will improve speed, but some matching sites near stronger ones (but beyond 
#' minDist) might be lost when more than two sites overlap. Not recommended.
#' @param BP Pass `BiocParallel::MulticoreParam(ncores, progressbar=TRUE)` to enable 
#' multithreading.
#' @param verbose Logical; whether to print additional progress messages (default on if 
#' not multithreading)
#'
#' @return A GRanges of all matches
#' 
#' @import BiocParallel Biostrings GenomicRanges
#' @export
#'
#' @examples
#' # we create mock RNA sequences and seeds:
#' seqs <- sapply(1:10, FUN=function(x) paste(sample(strsplit("ACGU", "")[[1]], 
#'                                      1000, replace=TRUE),collapse=""))
#' names(seqs) <- paste0("seq",1:length(seqs))
#' seeds <- c("AAACCAC", "AAACCUU")
#' findSeedMatches(seqs, seeds)
findSeedMatches <- function( seqs, seeds, seedtype=c("auto", "RNA","DNA"), shadow=0, 
                             minLogKd=1, types=NULL, max.noncanonical.motifs=Inf, 
                             keepMatchSeq=FALSE, minDist=7L, fastRemoveOverlaps=FALSE, 
                             BP=NULL, verbose=NULL){
  library(GenomicRanges)
  library(stringr)
  library(BiocParallel)
  library(Biostrings)
  if(is(seeds, "CompressedKdModelList")) seeds <- decompressKdModList(seeds)
  if(is.null(verbose)) verbose <- is(seeds,"KdModel") || length(seeds)==1 || is.null(BP)
  if(verbose) message("Preparing sequences...")
  args <- .prepSeqs(seqs, seeds, seedtype, shadow)
  seqs <- args$seqs
  if("seeds" %in% names(args)) seeds <- args$seeds
  rm(args)

  if(is(seeds,"KdModel") || length(seeds)==1){
    if(all(class(seeds)=="list")) seeds <- seeds[[1]]
    if(is.null(verbose)) verbose <- TRUE
    m <- .find1SeedMatches(seqs, seeds, keepMatchSeq, minDist=minDist, minLogKd=minLogKd,
                           types=types, max.noncanonical.motifs=max.noncanonical.motifs,
                           fastRemoveOverlaps=fastRemoveOverlaps, verbose=verbose)
  }else{
    if(is.null(BP)) BP <- SerialParam()
    if(is.null(verbose)) verbose <- !(bpnworkers(BP)>1 | length(seeds)>5)
    m <- bplapply( seeds, seqs=seqs, verbose=verbose, minDist=minDist, minLogKd=minLogKd,
                   types=types, max.noncanonical.motifs=max.noncanonical.motifs,
                   fastRemoveOverlaps=fastRemoveOverlaps,
                   FUN=.find1SeedMatches, BPPARAM=BP)
    m <- m[!sapply(m,is.null)]
    m <- unlist(GRangesList(m))
  }

  gc(verbose = FALSE, full = TRUE)
  
  row.names(m) <- NULL
  start(m) <- start(m)-3
  end(m) <- end(m)-3
  
  if(shadow>0){
    end(m) <- end(m)+shadow
    start(m) <- start(m)+shadow
  }
  names(m) <- NULL
  
  metadata(m)$shadow <- shadow
  metadata(m)$minDist <- minDist
  if(!is.null(types)) metadata(m)$types <- types
  metadata(m)$max.noncanonical.motifs <- max.noncanonical.motifs
  
  m
}


#' removeOverlappingMatches
#' 
#' Removes elements from a GRanges that overlap (or are within a given distance of) other 
#' elements higher up in the list (i.e. assumes that the ranges are sorted in order of
#' priority). The function handles overlaps between more than two ranges by successively
#' removing those that overlap higher-priority ones
#'
#' @param x A GRanges, sorted by (decreasing) importance
#' @param minDist Minimum distance between ranges
#'
#' @return A filtered GRanges.
#' @export
#'
#' @examples
removeOverlappingMatches <- function(x, minDist=7L, method=c("exact","fast")){
  method <- match.arg(method)
  red <- GenomicRanges::reduce(x, with.revmap=TRUE, min.gapwidth=minDist)$revmap
  red <- red[lengths(red)>1]
  if(length(red)==0) return(x)
  if(method=="fast"){
    toRemove <- setdiff(unlist(red),min(red))
  }else{
    red2 <- red[lengths(red)==2]
    toRemove <- setdiff(unlist(red2), min(red2))
    toRemove <- c(toRemove, unlist(lapply( red[lengths(red)>2], FUN=function(i){
      w <- j <- i <- sort(i)
      y <- ranges(x)[i]
      while(length(y)>1 && length(w)>0){
        # remove anything overlappnig the top one
        w <- 1+which(overlapsAny(y[-1],y,maxgap=minDist))
        if(length(w)>0){
          j <- j[-w]
          y <- y[-w]
        }
      }
      return(setdiff(i,j))
    })))
  }
  if(length(toRemove)>0) x <- x[-toRemove]
  x
}


.prepSeqs <- function(seqs, seeds, seedtype=c("auto", "RNA","DNA"), shadow=0){
  if(is.null(names(seqs))) names(seqs) <- paste0("seq",seq_along(seqs)) 
  seedtype <- match.arg(seedtype)
  seqtype <- .guessSeqType(seqs)
  ret <- list()
  if( is(seeds, "KdModel") || 
      (is.list(seeds) && all(sapply(seeds, FUN=function(x) is(x,"KdModel")))) ){
    if(is.null(names(seeds)))
      stop("If `seeds` is a list of kd models, it should be named.")
    if(seedtype=="RNA" || seqtype=="RNA") 
      stop("If `seeds` is a list of kd models, both the seeds and the target
sequences should be in DNA format.")
  }else{
    if(is.null(names(seeds))) n <- names(seeds) <- seeds
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
    ret$seeds <- seeds
  }
  seqnms <- names(seqs)
  if(is.character(seqs)) seqs <- DNAStringSet(seqs)
  names(seqs) <- seqnms
  seqs <- seqs[lengths(seqs)>=(shadow+7)]
  if(shadow>0) seqs <- subseq(seqs, shadow+1)
  seqs <- padAndClip(seqs, views=IRanges(start=1-3, width=lengths(seqs)+6), 
                     Lpadding.letter = "N", Rpadding.letter = "N")
  c(ret, list(seqs=seqs))
}

.find1SeedMatches <- function(seqs, seed, keepMatchSeq=FALSE, types=NULL, minLogKd=0,
                              max.noncanonical.motifs=Inf, minDist=1, 
                              fastRemoveOverlaps=FALSE, verbose=FALSE){
  library(GenomicRanges)
  library(stringr)
  library(BiocParallel)
  library(Biostrings)

  if(!is.null(types) && length(types)>0){
    types <- match.arg(types, c("8mer","7mer-m8","7mer-A1","6mer","offset 6mer",
                                "non-canonical"),
                       several.ok = TRUE)
  }else{
    types <- NULL
  }
    
  if(verbose) message("Scanning...")
  
  maxnc <- max.noncanonical.motifs
  if(!is.null(types) && !("non-canonical" %in% types)) maxnc <- 0
  m <- .scanOneSeed(seqs=seqs, seed=seed, max.noncanonical.motifs=maxnc)
  
  if(is.null(m)){
    if(verbose) message("Nothing found!")
    return(NULL)
  }
  
  m <- keepSeqlevels(m, seqlevelsInUse(m))
  m <- m[order(seqnames(m))]
  
  if(verbose) message("Extracting sequences...")
  
  seqs <- seqs[seqlevels(m)]

  r <- ranges(m)
  start(r) <- start(r)-3
  end(r) <- end(r)+3

  r <- split(r, seqnames(m))
  names(r) <- NULL
  
  ms <- unlist(extractAt(seqs, r))
  names(ms) <- NULL
  ms <- Rle(as.factor(ms))
  
  row.names(m) <- NULL
  
  if(verbose) message("Characterizing matches...")
  
  mc <- characterizeSeedMatches( ms, seed )
  mc$type <- Rle(mc$type)
  row.names(mc) <- NULL
  if(keepMatchSeq) m$sequence <- ms
  rm(ms)
  mcols(m) <- cbind(mcols(m), mc)
  rm(mc)
  if(!is.null(types)) m <- m[m$type %in% types,]
  if(!is.null(minLogKd) && minLogKd!=0 && "log_kd" %in% colnames(mcols(m))){
    if(minLogKd>0) minLogKd <- -minLogKd
    m <- m[m$log_kd <= minLogKd]
  }
  
  if(is(seed,"KdModel")){
    m <- m[order(seqnames(m), m$log_kd),]
    m <- m[order(m$log_kd),]
    m$seed <- Rle(rep(as.factor(seed$name),length(m)))
  }else{
    m <- m[order(seqnames(m), m$type),]
    m$seed <- Rle(rep(as.factor(seed),length(m)))
  }
  
  if(minDist>-Inf){
    if(verbose) message("Removing overlaps...")
    m <- removeOverlappingMatches(m, minDist=minDist, 
                                  method=ifelse(fastRemoveOverlaps,"fast","exact"))
  }
  
  names(m) <- NULL
  m
}


.scanOneSeed <- function(seqs, seed, max.noncanonical.motifs=Inf){
  seqnms <- names(seqs)
  seqs <- as.character(seqs)
  names(seqs) <- seqnms
  if(is(seed,"KdModel")){
    seed2 <- names(summary(seed))
    seedtypes <- sapply(seed2, seed=seed$canonical.seed, FUN=.getMatchType)
    if(max.noncanonical.motifs==0){
      seed2 <- seed2[seedtypes!="non-canonical"]
    }else{
      seed2.a <- seed2[seedtypes!="non-canonical"]
      seed2 <- seed2[seedtypes=="non-canonical"]
      seed2 <- seed2[seq_len(min(length(seed2), max.noncanonical.motifs))]
      seed2 <- c(seed2.a,seed2)
    }
    seed2 <- substring(seed2,2)
  }else{
    seed2 <- substr(seed,2,7)
  }
  # look-around matching to get overlapping seeds
  pos <- gregexpr( paste0("(?=",paste(unique(seed2),collapse="|"),")"),
                   seqs, perl=TRUE )
  pos <- lapply(lapply(pos, as.numeric), y=-1, setdiff)
  if(sum(sapply(pos,length))==0) return(NULL)
  GRanges( rep(names(seqs), sapply(pos,length)), 
           IRanges( start=unlist(pos), width=6 ) )
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
#' @param seed The seed which sequences of `x` should match, or a 'KdModel' object.
#'
#' @return A data.frame, or a `GRanges` if `x` is a `GRanges`
#' @export
#'
#' @examples
#' characterizeSeedMatches(c("UAAACCACCC","CGAACCACUG"), "AAACCAC")
characterizeSeedMatches <- function(x, seed=NULL){
  if(is(seed, "KdModel")){
    kd.model <- seed
    sseed <- seed$canonical.seed
  }else{
    kd.model <- NULL
    sseed <- seed
  }
  sseed <- as.character(sseed)
  if(length(sseed)>1 || nchar(sseed)!=7) 
    stop("`seed` must contain a single string of 7 characters.")
  if(!is.character(x) && !is.factor(x) && !is(x,"Rle")){
    y <- characterizeSeedMatches(x$sequence, seed)
    if(is(x,"GRanges")){
      mcols(x) <- cbind(mcols(x), y)
    }else{
      x <- cbind(x,y)
    }
    return(x)
  }
  if(is(x,"Rle") || any(duplicated(x))){
    # we resolve type for unique sequences
    if(!is.null(levels(x))){
      levels(x) <- gsub("N","x",levels(x))
      y <- characterizeSeedMatches(levels(x), seed)
    }else{
      x <- gsub("N","x",x)
      y <- characterizeSeedMatches(unique(x), seed)
    }
    y <- y[as.character(x),,drop=FALSE]
    row.names(y) <- NULL
    return(y)
  }
  d <- data.frame( row.names=x, 
                   type=vapply( as.character(x), seed=sseed, 
                                FUN=.getMatchType, FUN.VALUE=character(1) ) )
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

#' aggregateMatches
#'
#' @param e A GRanges object as produced by `findSeedMatches`.
#'
#' @return An aggregated data.frame
#' @import data.table
#' @importFrom GenomicRanges mcols
#' @export
aggregateMatches <- function(e){
  d <- as.data.frame(mcols(e))
  d$transcript <- as.factor(seqnames(e))
  d <- as.data.table(d)
  d2 <- subset(d, type!="non-canonical")
  m <- dcast( d2[,.(N=.N), by=c("transcript","seed","type")],
                 formula=transcript+seed~type, value.var="N", fill=0)
  ag2 <- ag1 <- NULL
  if("log_kd" %in% colnames(d)){
    ag1 <- d2[,.( log_kd.canonical=-log10(sum(10^-log_kd))), by=c("transcript","seed")]
    ag2 <- d[,.( log_kd=-log10(sum(10^-log_kd))), by=c("transcript","seed")]
  }
  rm(d,d2)
  if(!is.null(ag2)){
    m <- merge(m, ag1, by=c("transcript","seed"), all=TRUE)
    m <- merge(m, ag2, by=c("transcript","seed"), all=TRUE)
    m$log_kd.canonical[is.na(m$log_kd.canonical)] <- 0
  }
  m <- as.data.frame(m)
  m$`offset 6mer` <- NULL
  ff <- c("8mer","7mer-m8","7mer-A1","6mer")
  for( f in ff) m[[f]][is.na(m[[f]])] <- 0
  colnames(m)[3:6] <- gsub("-","",paste0("n.",colnames(m)[3:6]))
  m
}

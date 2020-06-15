#' IndexedFst
#' 
#' Objects of the IndexedFst class enable fast named random access to FST files.
#' This is particularly appropriate for large data.frames which often need to 
#'  be accessed according to the (e.g. factor) value of a particular column.
#' 
#' @examples
#' # we first create and save an indexed FST file
#' tmp <- tempdir()
#' f <- system.file(tmp, "test")
#' d <- data.frame( category=sample(LETTERS[1:4], 10000, replace=TRUE),
#'                  var2=sample(LETTERS, 10000, replace=TRUE),
#'                  var3=runif(10000) )
#' format(object.size(d1),units="Kb")
#' saveIndexedFst(d, "category", f)
#' rm(d)
#' # we then load the index, and can use category names for random access:
#' d <- loadIndexedFst(f)
#' format(object.size(d),units="Kb")
#' nrow(d)
#' names(d)
#' d$A
#' @seealso \code{\link{saveIndexedFst}}, \code{\link{loadIndexedFst}}
#' @export
#' @import methods
#' @exportClass IndexedFst
#' @author Pierre-Luc Germain, \email{pierre-luc.germain@@hest.ethz.ch}
#' @aliases IndexedFst-class IndexedFst
#' @rdname IndexedFst-class
#' @name IndexedFst-class
setClass(
  "IndexedFst",
  slots=representation(
    fst.file="character",
    index="data.frame",
    nthreads="integer"
  ),
  prototype=prototype(fst.file=NA_character_, index=data.frame(), nthreads=1L),
  validity=function(object){
    if(length(object@nthreads)!=1 || !is.integer(object@nthreads) || object@nthreads<1)
      stop("`nthreads` should be a positive integer")
    if(length(object@fst.file)!=1) stop("fst.file should be a character of length 1.")
    if(!file.exists(object@fst.file)) stop("FST file does not exist!")
    TRUE
  }
)

setMethod("initialize", "IndexedFst", function(.Object, ...) {
  o <- callNextMethod(.Object, ...)
  o@fst.file <- normalizePath(o@fst.file)
  o@nthreads <- as.integer(o@nthreads)
  ff <- gsub("\\.fst$",".idx",o@fst.file)
  o@index <- tryCatch( read.delim(ff, header=FALSE, row.names=1),
                       error=function(e) stop("Could not find or read index file."))
  validObject(o)
  return(o)
})


#' @rdname IndexedFst-class
#' @importMethodsFrom methods show
#' @export
setMethod("show", "IndexedFst", function(object){
  paste0(ff@fst.file, " (",nrow(ff@index)," sets)")
})

#' @rdname IndexedFst-class
#' @importFrom fst fst.metadata
#' @export
setMethod("summary", "IndexedFst", function(object){
  fst.metadata(object@fst.file)
})

#' @rdname IndexedFst-class
#' @export
setMethod("names", "IndexedFst", function(x){
  row.names(x@index)
})

#' @rdname IndexedFst-class
#' @export
setMethod("length", "IndexedFst", function(x){
  nrow(x@index)
})

#' @rdname IndexedFst-class
#' @export
setMethod("lengths", "IndexedFst", function(x){
  y <- x@index[,2]-x@index[,1]+1
  names(y) <- names(x)
  y
})

setGeneric("nrow")
#' @rdname IndexedFst-class
#' @export
setMethod("nrow", "IndexedFst", function(x){
  max(x@index[,2])
})

setGeneric("ncol")
#' @rdname IndexedFst-class
#' @export
setMethod("ncol", "IndexedFst", function(x){
  length(fst.metadata(x@fst.file)$columnNames)
})

setGeneric("colnames")
#' @rdname IndexedFst-class
#' @export
setMethod("colnames", signature("IndexedFst"), function(x){
  fst.metadata(x@fst.file)$columnNames
})

#' @rdname IndexedFst-class
#' @export
setMethod("[[", signature("IndexedFst"), function(x, i, j=NULL, ...){
  if(is.numeric(i)){
    name <- names(x)[i]
  }else{
    name <- i
  }
  .fst.read.wrapper(x, name)
})

#' @rdname IndexedFst-class
#' @export
setMethod("[", signature("IndexedFst"), function(x, i, j=NULL, ...){
  if(is.logical(i)) i <- which(i)
  if(is.numeric(i)){
    name <- names(x)[i]
  }else{
    name <- i
  }
  .fst.read.wrapper(x, name)
})

#' @rdname IndexedFst-class
#' @export
setMethod("$", "IndexedFst", definition = function(x, name){
  .fst.read.wrapper(x, match.arg(name, row.names(x@index)))
})

#' @rdname IndexedFst-class
#' @export
setMethod("head", "IndexedFst", definition = function(x, n=6L, ...){
  if(!is.numeric(n) || !(n>0) || n!=as.integer(n))
    stop("`n` should be a positive integer.")
  w <- which(x@index[,2]>=n)
  if(length(w)==0) return(.fst.read.wrapper(x))
  head(x[seq_len(w[1])], n=n)
})

#' @rdname IndexedFst-class
#' @export
setMethod("as.data.frame", "IndexedFst", definition=function(x, name) {
  .fst.read.wrapper(x)
})

#' @importFrom fst threads_fst read.fst
#' @importFrom data.table rbindlist 
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
.fst.read.wrapper <- function(x, names=NULL, convertGR=TRUE){
  if(!is(x,"IndexedFst")) stop("`x` should be an IndexedFst object.")
  if(!is.null(names) && length(names)>1){
    f <- lapply(names, FUN=function(i) .fst.read.wrapper(x,i,convertGR=FALSE))
    if(length(f)>5){
      f <- as.data.frame(data.table::rbindlist(f))
    }else{
      f <- do.call(rbind, f)
    }
  }else{
    ont <- threads_fst()
    threads_fst(x@nthreads)
    if(is.null(names)){
      f <- read.fst(x@fst.file)
    }else{
      f <- read.fst(x@fst.file, from=x@index[names,1], to=x@index[names,2])
    }
    threads_fst(ont)
  }
  if(convertGR){
    if(all(c("seqnames","start","end") %in% colnames(f))){
      f <- GRanges(seqnames=f$seqnames, IRanges(f$start, f$end), strand=f$strand,
                   f[,setdiff(colnames(f), c("seqnames","start","end","width","strand"))])
    }
  }
  f
}

#' loadIndexedFst
#' 
#' Loads the index of an FST file, for fast random access.
#'
#' @param file Path to the fst file, it's index (.idx), or their prefix.
#' @param nthreads Number of threads to use for reading (default 1). This does not affect
#' the loading of the index itself, but will affect all downstream reading operations
#' performed on the object.
#'
#' @return An object of class IndexedFst
#' @seealso \code{\link{IndexedFst-class}}
#' @export
loadIndexedFst <- function(file, nthreads=1){
  if(grepl("\\.fst$",file)) return(new("IndexedFst", fst.file=file))
  if(grepl("\\.idx$",file)) return(new("IndexedFst", fst.file=gsub("idx$","fst",file)))
  new("IndexedFst", fst.file=paste0(file,".fst"), nthreads=as.integer(nthreads))
}

#' saveIndexedFst
#' 
#' Saves a data.frame into an indexed FST file.
#'
#' @param d A data.frame
#' @param index.by A column of `d` by which it should be indexed.
#' @param file.prefix Path and prefix of the output files.
#' @param nthreads Number of threads to use when writing (default 1). If NULL, will use
#' `threads_fst()`
#' @param ... Passed to `write.fst`
#'
#' @return Nothing.
#' @seealso \code{\link{IndexedFst-class}}
#' @importFrom fst write.fst threads_fst
#' @export
saveIndexedFst <- function(d, index.by, file.prefix, nthreads=1, ...){
  if(!is.data.frame(d)) stop("`d` should be a data.frame.")
  if((!is.character(index.by) && !is.integer(index.by)) || length(index.by)!=1 ||
     is.null(d[[index.by]]))  stop("`index.by` should be a scalar character or integer ",
                                   "indicating the column by which to index.")
  if(is.numeric(d[[index.by]]))
    warning("The `index.by` column selected is numeric! ",
            "If the column is really continuous, rather than having many repeated ",
            "values, this might lead to very suboptimal behaviors!")
  if(!is.null(nthreads)) threads_fst(nthreads)
  d <- d[order(d[[index.by]]),]
  file.prefix <- gsub("\\.fst$","",file.prefix)
  write.fst(d, paste0(file.prefix, ".fst"), ...)
  idx <- lapply(split(seq_len(nrow(d)), d[[index.by]]), FUN=range)
  idx <- data.frame( seed=names(idx), start=sapply(idx,FUN=function(x) x[1]), 
                      end=sapply(idx,FUN=function(x) x[2]))
  write.table(idx, paste0(file.prefix, ".idx"), row.names=FALSE, col.names=FALSE, 
              sep="\t",quote=F)
}
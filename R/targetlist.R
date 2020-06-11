#' @export
#' @importFrom fst read.fst
setClass(
  "IndexedFst",
  slots=representation(
    fst.file="character",
    index="data.frame"
  ),
  prototype=prototype(fst.file=NA_character_, index=data.frame()),
  validity=function(object){
    if(length(object@fst.file)!=1) stop("fst.file should be a character of length 1.")
    if(!file.exists(object@fst.file)) stop("FST file does not exist!")
    TRUE
  }
)

setMethod("initialize", "IndexedFst", function(.Object, ...) {
  o <- callNextMethod(.Object, ...)
  o@fst.file <- normalizePath(o@fst.file)
  ff <- gsub("\\.fst$",".idx",o@fst.file)
  o@index <- tryCatch( read.delim(ff, header=FALSE, row.names=1),
                       error=function(e) stop("Could not find or read index file."))
  validObject(o)
  return(o)
})


#' @export
setMethod("show", "IndexedFst", function(object){
  paste0(ff@fst.file, " (",nrow(ff@index)," sets)")
})

#' @export
setMethod("names", "IndexedFst", function(x){
  row.names(x@index)
})

#' @export
setMethod("length", "IndexedFst", function(x){
  nrow(x@index)
})

#' @export
setMethod("lengths", "IndexedFst", function(x){
  y <- x@index[,2]-x@index[,1]+1
  names(y) <- names(x)
  y
})

#' @export
setMethod("nrow", "IndexedFst", function(x){
  max(x@index[,2])
})

#' @export
setMethod("[[", signature(x = "IndexedFst"), function(x, i, j=NULL, ...){
  if(is.numeric(i)){
    name <- names(x)[i]
  }else{
    name <- i
  }
  read.fst(x@fst.file, from=x@index[i,1], to=x@index[i,2])
})

#' @export
setMethod("[", signature(x = "IndexedFst"), function(x, i, j=NULL, ...){
  if(is.logical(i)) i <- which(i)
  if(is.numeric(i)){
    name <- names(x)[i]
  }else{
    name <- i
  }
  do.call(rbind, lapply(name, FUN=function(i) x[[i]]))
})

#' @export
setMethod("$", "IndexedFst", definition = function(x, name) {
  name <- match.arg(name, row.names(x@index))
  read.fst(x@fst.file, from=x@index[name,1], to=x@index[name,2])
})

#' @export
setMethod("as.data.frame", "IndexedFst", definition=function(x, name) {
  fst.read(x@fst.file)
})


loadIndexedFst <- function(file){
  if(grepl("\\.fst$",file)) return(new("IndexedFst", fst.file=file))
  if(grepl("\\.idx$",file)) return(new("IndexedFst", fst.file=gsub("idx$","fst",file)))
  new("IndexedFst", fst.file=paste0(file,".fst"))
}

saveIndexedFst <- function(d, index.by, file.prefix, ...){
  if(!is.data.frame(d)) stop("`d` should be a data.frame.")
  if((!is.character(index.by) && !is.integer(index.by)) || length(index.by)!=1 ||
     is.null(d[[index.by]]))  stop("`index.by` should be a scalar character or integer ",
                                   "indicating the column by which to index.")
  d <- d[order(d[[index.by]]),]
  file.prefix <- gsub("\\.fst$","",file.prefix)
  write.fst(m, paste0(file.prefix, ".fst"), ...)
  idx <- lapply(split(seq_len(nrow(d)), d[[index.by]]), FUN=range)
  idx <- data.frame( seed=names(idx), start=sapply(idx,FUN=function(x) x[1]), 
                      end=sapply(idx,FUN=function(x) x[2]))
  write.table(idx, paste0(file.prefix, ".idx"), row.names=FALSE, col.names=FALSE, sep="\t",quote=F)
}
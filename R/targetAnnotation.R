setClass(
  "targetAnnotation",
  slots=representation(
    object="list",
    set.properties="data.frame",
    metadata="list",
    feature.synonyms="factor",
    sets.synonyms="factor"
  ),
  validity=function(object){
    stopifnot(length(object@object)==1)
    TRUE
  }
)

#' @rdname targetAnnotation
#' @importMethodsFrom methods show
#' @export
setMethod("show", "targetAnnotation", function(object){
  paste0("A target annotation (",class(object),") containing ", length(object),
         " sets of among ", length(setTargetFeatures(object)), " targets")
})

#' @rdname targetAnnotation
#' @export
setClass(
  "targetAnnotationMatrix",
  contains="targetAnnotation",
  slots=c(content="character"),
  prototype=prototype(metadata=list(), content="score"),
  validity=function(object){
    stopifnot(length(object@object)==1)
    stopifnot(length(object@content)==1)
    stopifnot(.is.matrix(object@object[[1]]))
    stopifnot(is.logical(object@object[[1]][1,1]) || 
                is.numeric(object@object[[1]][1,1]))
    stopifnot(all(row.names(object@set.properties)==colnames(object@object[[1]])))
    stopifnot(!is.null(row.names(object@object[[1]])) &&
                !is.null(colnames(object@object[[1]])))
    stopifnot(all(levels(object@feature.synonyms) %in% 
                    row.names(object@object[[1]])))
    stopifnot(all(levels(object@sets.synonyms) %in% 
                    colnames(object@object[[1]])))
    TRUE
  }
)

#' @rdname targetAnnotation
#' @export
setClass(
  "targetAnnotationDF",
  contains="targetAnnotation",
  validity=function(object){
    stopifnot(length(object@object)==1)
    stopifnot(is.data.frame(object@object[[1]]) || 
                is(object@object[[1]], "DFrame"))
    stopifnot(all(c("feature","set") %in% colnames(object@object[[1]])))
    stopifnot(all(row.names(object@set.properties)==levels(object@object[[1]]$set)))
    stopifnot(is.factor(object@object[[1]]$feature))
    stopifnot(is.factor(object@object[[1]]$set))
    stopifnot(all(levels(object@feature.synonyms) %in% 
                    levels(object@object[[1]]$feature)))
    stopifnot(all(levels(object@sets.synonyms) %in% 
                    levels(object@object[[1]]$set)))
    TRUE
  }
)


#' @rdname targetAnnotation
#' @importMethodsFrom methods show
#' @export
setMethod("show", "targetAnnotationDF", function(object){
  paste0("A target annotation (",class(object),") containing ", length(object),
         " sets of among ", length(setTargetFeatures(object)), " targets, ",
         "with the attributes:\n", 
         paste(setdiff(colnames(object@object[[1]]), c("set","feature")),
               collapse=", "))
})
#' @rdname targetAnnotation
#' @export
setMethod("head", "targetAnnotationDF", definition = function(x, n=6L, ...){
  head(x@object[[1]], n)
})


#' @rdname targetAnnotation
#' @export
setMethod("as.data.frame", "targetAnnotationDF", definition=function(x, name){
  x@object[[1]]
})
#' @rdname targetAnnotation
#' @export
setMethod("as.data.frame", "targetAnnotationMatrix", definition=function(x, name){
  content <- x@content
  x <- x@object[[1]]
  if(is(x,"sparseMatrix")){
    dimn <- dimnames(x)
    x <- summary(x)
    w <- which(x$x!=ifelse(is.logical(x$x),FALSE,0))
    xout <- data.frame( feature=factor(x$i[w], levels=seq_len(length(dimn[[1]])), labels=dimn[[1]]),
                        set=factor(x$j[w], levels=seq_len(length(dimn[[2]])), labels=dimn[[2]]) )
    if(!is.logical(x$x)) xout[[content]] <- x$x[w]
  }else{
    w <- which(as.vector(x!=ifelse(is.logical(x[1,1]),FALSE,0)))
    xout <- data.frame( feature=factor(rep(seq_len(nrow(x)),ncol(x))[w], levels=seq_len(nrow(x)), labels=row.names(x)),
                        set=factor(rep(seq_len(ncol(x)),each=nrow(x))[w], levels=seq_len(ncol(x)), labels=row.names(x)))
    if(!is.logical(x[1,1])) xout[[content]] <- as.vector(x)[w]
  }
  xout
})

#' @rdname targetAnnotation
#' @export
setMethod("as.matrix", "targetAnnotationMatrix", definition=function(x, ...){
  x@object[[1]]
})

#' @rdname targetAnnotation
#' @export
setMethod("as.matrix", "targetAnnotationDF", definition=function(x, ...){
  x <- x@object[[1]]
  feats <- row.names(x)
  col <- head(intersect(c("score","sites"), colnames(x)),1)
  .setsToScoreMatrix(setNames(seq_along(feats), feats), 
                     sets=x, column=col, keepSparse=TRUE)
})

.ob <- function(x) x@object[[1]]

#' @export
setGeneric("setTargetFeatures", function(x){ })
#' @rdname targetAnnotation
#' @export
setMethod("setTargetFeatures", "targetAnnotationDF", function(x){
  levels(x@object[[1]]$features)
})
#' @rdname targetAnnotation
#' @export
setMethod("setTargetFeatures", "targetAnnotationMatrix", function(x){
  row.names(x@object[[1]])
})

#' @rdname targetAnnotation
#' @export
setMethod("names", "targetAnnotationDF", function(x){
  levels(x@object[[1]]$set)
})
#' @rdname targetAnnotation
#' @export
setMethod("names", "targetAnnotationMatrix", function(x){
  colnames(x@object[[1]])
})

#' @rdname targetAnnotation
#' @export
setMethod("length", "targetAnnotationMatrix", function(x){
  ncol(x@object[[1]])
})
#' @rdname targetAnnotation
#' @export
setMethod("length", "targetAnnotationDF", function(x){
  length(levels(x@object[[1]]$set))
})



#' @rdname targetAnnotation
#' @export
setMethod("[", signature("targetAnnotationDF"), function(x, i, j=NULL, ...){
  if(!missing(i)){
    if(is.logical(i)) i <- which(i)
    if(is.factor(i)) i <- as.character(i)
    if(is.character(i)) i <- as.integer(factor(i,setTargetFeatures(x)))
    i <- as.integer(i)
    i <- unique(i[!is.na(i)])
    x@feature.synonyms <- droplevels(x@feature.synonyms[x@feature.synonyms %in% names(x)[i]])
    x@object[[1]] <- x@object[[1]][as.integer(x@object[[1]]$feature) %in% i,,drop=FALSE]
    x@object[[1]]$feature <- factor(x@object[[1]]$feature,setTargetFeatures(x)[i])
  }
  if(!missing(j)){
    if(is.logical(j)) j <- which(j)
    if(is.factor(j)) j <- as.character(j)
    if(is.character(j)) j <- as.integer(factor(j,names(x)))
    j <- as.integer(j)
    j <- unique(j[!is.na(j)])
    x@set.properties <- x@set.properties[j,]
    x@sets.synonyms <- droplevels(x@sets.synonyms[x@sets.synonyms %in% names(x)[j]])
    x@object[[1]] <- x@object[[1]][as.integer(x@object[[1]]$set) %in% j,,drop=FALSE]
    x@object[[1]]$set <- factor(x@object[[1]]$set,names(x)[j])
  }
  x
})

#' @rdname targetAnnotation
#' @export
setMethod("[", signature("targetAnnotationMatrix"), function(x, i, j=NULL, ...){
  if(!missing(i)){
    if(is.logical(i)) i <- which(i)
    if(is.factor(i)) i <- as.character(i)
    if(is.character(i)) i <- as.integer(factor(i,setTargetFeatures(x)))
    if(is.numeric(i)) i <- as.integer(i)
    i <- unique(i[!is.na(i)])
    x@feature.synonyms <- droplevels(x@feature.synonyms[x@feature.synonyms %in% names(x)[i]])
    x@object[[1]] <- x@object[[1]][i,,drop=FALSE]
  }
  if(!missing(j)){
    if(is.logical(j)) j <- which(j)
    if(is.factor(j)) j <- as.character(j)
    if(is.character(j)) j <- as.integer(factor(j,names(x)))
    if(is.numeric(j)) j <- as.integer(j)
    j <- unique(j[!is.na(j)])
    x@set.properties <- x@set.properties[j,]
    x@sets.synonyms <- droplevels(x@sets.synonyms[x@sets.synonyms %in% names(x)[j]])
    x@object[[1]] <- x@object[[1]][,j,drop=FALSE]
  }
  x
})



#' targetAnnotation
#' 
#' Creates a new targetAnnotation object.
#'
#' @param x A names list of items, or a (sparse) matrix with named dimensions,
#' or a data.frame with at least the columns `feature` and `set`, or a 
#' (h5/rds) file containing one of these.
#' @param content The content of the matrix (e.g. 'score' or 'sites'); ignored
#' if `x` is a data.frame or list.
#' @param metadata A list of metadata.
#' @param set.properties A data.frame of additional set properties.
#' @param feature.synonyms A named factor vector of feature synonyms.
#' @param sets.synonyms A named factor vector of set names synonyms.
#'
#' @return An object of class `targetAnnotation`; depending on `x`, either a 
#' `targetAnnotationMatrix` or a `targetAnnotationDF`.
#' @rdname targetAnnotation
#' @importFrom HDF5Array HDF5Array
#' @import DelayedArray
#' @export
targetAnnotation <- function(x, content="score", metadata=list(), set.properties=NULL,
                             feature.synonyms=NULL, sets.synonyms=NULL){
  if(is.character(x)){
    stopifnot(length(x)==1)
    if(grepl("\\.h5$", x)){
      x <- HDF5Array(x,as.sparse=TRUE)
    }else if(grepl("\\.rds$", x)){
      x <- readRDS(x)
    }else{
      stop("File format not yet implemented, please read in manually first.")
    }
  }
  if(.is.matrix(x)){
    if(is.null(set.properties)) set.properties <- data.frame(row.names=colnames(x))
    if(is.null(feature.synonyms)) feature.synonyms <- factor(c(), levels=row.names(x))
    if(is.null(sets.synonyms)) sets.synonyms <- factor(c(), levels=colnames(x))
    x <- new("targetAnnotationMatrix", object=list(x), content=content, 
             metadata=metadata, set.properties=set.properties, 
             feature.synonyms=feature.synonyms, sets.synonyms=sets.synonyms)
  }else{
    x <- .list2DF(x)
    stopifnot(all(c("feature","set") %in% colnames(x)))
    x$feature <- as.factor(x$feature)
    x$set <- as.factor(x$set)
    if(is.null(set.properties)) set.properties <- data.frame(row.names=levels(x$set))
    if(is.null(feature.synonyms)) feature.synonyms <- factor(c(), levels=levels(x$feature))
    if(is.null(sets.synonyms)) sets.synonyms <- factor(c(), levels=levels(x$set))
    x <- new("targetAnnotationDF", object=list(x), metadata=metadata, 
             set.properties=set.properties, feature.synonyms=feature.synonyms, 
             sets.synonyms=sets.synonyms)
  }
  x
}


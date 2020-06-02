#' @export
setClass(
  "CompressedKdModelList",
  slots=representation(
    frame="lm",
    sr="data.frame",
    fl="matrix",
    permir="list",
    classVersion="integer",
    metadata="list"
  ),
  prototype=prototype(frame=model.frame(~1), sr=data.frame(), fl=matrix(), 
                      permir=list(), classVersion=1L, 
                      metadata=list(created=Sys.Date())),
  validity=function(object){
    return(TRUE)
  }
)


#' @export
setMethod("as.list", "CompressedKdModelList", function (x){
  decompressKdModList(x)
})

#' @export
setMethod("show", "CompressedKdModelList", function(object){
  message("A `CompressedKdModelList` object (v", object@classVersion, 
  ") containing\nbinding affinity models from ", length(object@permir), 
  " miRNAs.")
})

#' @export
setMethod("summary", "CompressedKdModelList", function(object){
  show(object)
  if(!is.null(object@metadata$description)) message(object@metadata$description)
  if(!is.null(object@metadata$created)) message(object@metadata$created)
})

#' @export
setMethod("names", "CompressedKdModelList", function(x){
  names(x@permir)
})

#' @export
setMethod("metadata", "CompressedKdModelList", function (x) x@metadata )
#' @export
setMethod("metadata<-", "CompressedKdModelList", function (x, value){
  x@metadata <- value
  x
})

#' @export
setMethod("[[", signature(x = "CompressedKdModelList"), function(x, i, j=NULL, ...){
  name <- names(x)[i]
  .inflateKdMod(x@frame, co=x@sr[x@sr$miRNA==name,],
                o=x@permir[[name]],  x@fl[,name,drop=FALSE])
})

#' @export
setMethod("$", "CompressedKdModelList", definition = function(x, name) {
  name <- match.arg(name, names(x@permir))
  .inflateKdMod(x@frame, co=x@sr[x@sr$miRNA==name,],
                o=x@permir[[name]],  x@fl[,name,drop=FALSE])
})

#' decompressKdModList
#' 
#' @param mods A `CompressedKdModelList`
#'
#' @return A list of `KdModels`
#' @export
decompressKdModList <- function(mods){
  if(length(mods)==1 && is.character(mods)) mods <- readRDS(mods)
  if(!is(mods,"CompressedKdModelList")) stop("`mods` is not a CompressedKdModelList")
  SR <- split(mods$sr, mods$sr$miRNA)
  names(nn) <- nn <- names(mods$other)
  lapply(nn, FUN=function(n){
    .inflateKdMod(mods@frame, co=mods@sr[mods@sr$miRNA==n,],
                  o=mods@permir[[n]],  mods@fl[,n,drop=FALSE])
  })
}



#' compressKdModList
#'
#' @param mods A named list of objects of class `KdModel`
#'
#' @return A `CompressedKdModelList`
#' @export
compressKdModList <- function(mods){
  if(length(mods)==1 && is.character(mods)) mods <- readRDS(mods)
  if(!is.list(mods) || !all(sapply(mods, class2="KdModel", FUN=is)))
    stop("mods should be a named list of 'KdModel's!") 
  fl <- sapply(mods, FUN=function(x){
    co <- x$coefficients
    co <- co[grep("^fl",names(co))]
    x <- as.integer(round(100*co))
    names(x) <- names(co)
    x
  })
  sr <- dplyr::bind_rows(lapply(mods, FUN=function(x){
    co <- x$coefficients
    co <- co[c(1,grep("^sr|^ATRUE",names(co)))]
    col <- length(co)/2
    if(!all( names(co)[seq.int(from=2,to=col)] ==
             gsub(":ATRUE","",names(co)[seq.int(from=col+2,to=length(co))],fixed=TRUE)))
      stop("Model's coefficients are not standard!")
    data.frame( seed=names(co)[seq_len(col)],
                seed.coef=as.integer(round(100*co[seq_len(col)])),
                seed.A=as.integer(round(100*co[seq.int(col+1,length(co))])),
                stringsAsFactors = FALSE
    )
  }), .id="miRNA")
  sr$seed <- factor(sr$seed)
  sr$miRNA <- factor(sr$miRNA)
  
  names(otherfields) <- otherfields <- .kdmodfields()
  other <- lapply(mods, FUN=function(x){
    y <- lapply(otherfields, FUN=function(f) x[[f]])
    if(identical(y$qr$pivot, seq_len(length(y$qr$pivot))))
      y$qr$pivot <- length(y$qr$pivot)
    y$srlvls <- x$xlevels$sr
    y
  })
  
  modf <- mods[[1]]
  modf$xlevels <- list(sr=c(), fl=modf$xlevels$fl)
  for(f in otherfields) modf[[f]] <- NULL
  class(modf) <- "lm"
  
  new("CompressedKdModelList", frame=modf, sr=sr, fl=fl, permir=other)
}


.kdmodfields <- function(){
  c( "rank","qr","df.residual","mirseq","canonical.seed","pwm",
     "cor.with.cnn","mae.with.cnn","name" )
}
.inflateKdMod <- function(mod, co, o, fl){
  for(f in .kdmodfields()) if(f %in% names(o)) mod[[f]] <- o[[f]]
  if(length(mod$qr$pivot)==1) mod$qr$pivot <- seq_len(mod$qr$pivot)
  mod$xlevels$sr <- o$srlvls
  co2 <- c(co[,3],co[,4])/100
  names(co2) <- c(as.character(co$seed),paste0(co$seed,":ATRUE"))
  names(co2)[nrow(co)+1] <- "ATRUE"
  fl2 <- fl[,1]/100
  names(fl2) <- row.names(fl)
  mod$coefficients <- c( co2[grep(":",names(co2),invert=TRUE)],
                         fl2, co2[grep(":",names(co2))] )
  mod
}

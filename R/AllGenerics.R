#' @export
setClass(
	"enrich.results",
	slots=representation(
      input="list",
  		binary.signatures="list",
  		overlaps="list",
  		res="list",
  		info="list",
  		created="Date"
		),
    prototype=prototype( res=list(), info=list(objvers=2), input=list(), 
                         binary.signatures=list(), overlaps=list(), 
                         created=Sys.Date() )
)

setMethod("initialize", "enrich.results", function(.Object, ...) {
    o <- callNextMethod(.Object, ...)
    validObject(o)
    return(o)
})

setMethod("as.data.frame", "enrich.results", function (x, row.names = NULL, optional = FALSE, ...){
  x <- getResults(x)
  if(!is.null(row.names)) row.names(x) <- row.names
  x
})

#' @export
setMethod("show", "enrich.results", function(object){
    message("An `enrich.results` object with the following analyses:")
    print(lapply(object@res, FUN=function(x){
        head(x[,grep("set|pvalue|FDR|enrichment",colnames(x)),drop=FALSE])
    }))
})

#' @export
setMethod("summary", "enrich.results", function(object){
  message(paste("An `enrich.results` object with",length(object@res),"tests.\nAggregated top results:"))
  res <- getResults(object)
  if(!("enrichment" %in% colnames(res))){
    res$enrichment <- apply(abs(log2(as.matrix(res[,grep("enrichment",colnames(res)),drop=F])+0.05)),1,FUN=median)
  }
  res$medianP <- apply(res[,grep("pvalue",colnames(res)),drop=F],1,na.rm=TRUE, FUN=median)
  fields <- intersect(c("set","overlap","enrichment","medianP"),colnames(res))
  if(length(f <- grep("FDR", colnames(res), value=TRUE))>0)
    fields <- c(fields, f[1])
  head(res[order(res$medianP),fields])
})

#' @export
setMethod("$", "enrich.results", definition = function(x, name) {
  getResults(x,name)
})

#' @export
setMethod("names", "enrich.results", definition = function(x) {
  names(x@res)
})
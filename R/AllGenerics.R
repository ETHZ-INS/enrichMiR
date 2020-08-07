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
        head(x[,grep("miRNA|pvalue|FDR|enrichment",colnames(x)),drop=FALSE])
    }))
})

#' @export
setMethod("summary", "enrich.results", function(object){
  message(paste("An `enrich.results` object with",length(object@res),"tests.\nAggregated top results:"))
  res <- getResults(object)
  if(!("enrichment" %in% colnames(en))){
    if("EN.combined.enrichment" %in% colnames(res)){
      res$enrichment <- res$EN.combined.enrichment
      res$overlap <- res$EN.combined.overlap
      min.enr.thres <- -Inf
    }else{
      res$enrichment <- apply(abs(log2(as.matrix(res[,grep("enrichment"),drop=F])+0.05)),1,FUN=median)
    }
  }
  res$medianP <- apply(res[,grep("pvalue",colnames(res)),drop=F],1,FUN=median)
  head(res[order(res$medianP),intersect(c("miRNAs","expression","overlap","enrichment","medianP"),colnames(res))])
})

#' @export
setMethod("$", "enrich.results", definition = function(x, name) {
  getResults(x,name)
  }
)
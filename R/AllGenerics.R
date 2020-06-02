#' @export
setClass(
	"enrichMiR",
	slots=representation(
    DEA="data.frame",
		TS="data.frame",
		res="list",
		info="list",
		miRNA.expression="list",
		families="character",
		created="Date"
		),
    prototype=prototype(res=list(), info=list(objvers=1), miRNA.expression=list(family=NULL, miRNA=NULL), families=NA_character_, created=Sys.Date()),
    validity=function(object){
        if(!all(c("logFC","FDR") %in% colnames(object@DEA))) stop("The `DEA` slot should have at least the following columns: 'logFC', 'FDR'.")
        if(!all(c("family","feature") %in% colnames(object@TS))) stop("The `TS` slot should have at least the following columns: 'family', 'feature'.")
        #if(!any(as.character(object@TS$feature) %in% row.names(object@DEA))) stop("The row names of the `DEA` slot don't seem to match the `feature` column of the `TS` slot!")
        return(TRUE)
    }
)

setMethod("initialize", "enrichMiR", function(.Object, ...) {
    o <- callNextMethod(.Object, ...)
    validObject(o)
    return(o)
})

setMethod("as.data.frame", "enrichMiR", function (x, row.names = NULL, optional = FALSE, ...){
  x <- enrichMiR.results(x)
  if(!is.null(row.names)) row.names(x) <- row.names
  x
})

#' @export
setMethod("show", "enrichMiR", function(object){
    message("An `enrichMiR` object with the following analyses:")
    print(lapply(object@res, FUN=function(x){
        head(x[,grep("miRNA|pvalue|FDR|enrichment",colnames(x)),drop=FALSE])
    }))
})

#' @export
setMethod("summary", "enrichMiR", function(object){
  message(paste("An `enrichMiR` object with",length(object@res),"tests.\nAggregated top results:"))
  res <- enrichMiR.results(object)
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
setMethod("$", "enrichMiR", definition = function(x, name) {
    enrichMiR.results(x,name)
  }
)
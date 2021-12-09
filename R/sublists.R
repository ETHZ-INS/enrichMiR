#' getGOgenes
#'
#' @param go_ids A vector of GO IDs.
#' @param species A two-letter code indicating the species, such as `Mm` (default) or `Hs`.
#' @param translate.ids Logical; whether to translate GO IDs into terms (default TRUE)
#' @param ensembl_ids Logical; whether to translate Gene Symbols into Ensembl IDs
#'
#' @return
#' @export
#' @import AnnotationDbi GO.db
# @importFrom AnnotationDbi mget Term
getGOgenes <- function(go_ids, species="Mm", translate.ids=TRUE,ensembl_ids = FALSE){
  library(AnnotationDbi)
  if(!all(grepl("^GO:",go_ids)) && !is.null(names(go_ids))) go_ids <- names(go_ids)
  db <- paste0('org.',species,'.eg')
  library(package=paste0(db,'.db'), character.only = T)
  env <- get(paste0(db,'GO2ALLEGS'))
  eg <- mget(as.character(go_ids), env, ifnotfound=NA)
  x <- lapply(eg, FUN=function(x){ 
    x <- x[which(!is.na(x))]
    if(length(x)==0) return(c())
    if(ensembl_ids){ n <- 'ENSEMBL' }else{ n <- 'SYMBOL' }
    unique(as.character(unlist(mget(as.character(x), get(paste0(db,n))))))
  })
  x <- x[which(sapply(x,length)>0)]
  if(translate.ids) names(x) <- Term(names(x))
  if(ensembl_ids){x <- lapply(x, function(y) y[!is.na(y)])}
  x
}


#' findGO
#'
#' Finds GO categories where the term matches a query expression.
#'
#' @param expr The expression to query, interpreted as a regular expression if 
#' `fixed=FALSE`.
#' @param fixed Whether to interpret `expr` as a fixed expression, rather than 
#' a regular expression (default FALSE)
#' @param ontology Ontology to search, default all 3 (`BP|MF|CC`)
#'
#' @return A named vector of GO IDs, with Terms as names.
#' @export
# @import GO.db
findGO <- function(expr, fixed=FALSE, ontology="BP|MF|CC"){
  library("GO.db")
  a <- as.data.frame(GO.db::GOTERM)
  w <- intersect(grep(expr, a$Term, fixed=fixed, ignore.case=T), 
                 grep(ontology,a$Ontology))
  if(length(w)==0){
    warning("Nothing found!")
    return(NULL)
  }
  x <- a[w,c(1,3)]
  x <- x[!duplicated(x),]
  y <- x[,2]
  names(y) <- x[,1]
  return(y)
}

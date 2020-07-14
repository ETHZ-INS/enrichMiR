#' getGOgenes
#'
#' @param go_ids A vector of GO IDs.
#' @param species A two-letter code indicating the species, such as `Mm` (default) or `Hs`.
#' @param translate.ids Logical; whether to translate GO IDs into terms (default TRUE)
#'
#' @return
#' @export
getGOgenes <- function(go_ids, species="Mm", translate.ids=TRUE){
  library(AnnotationDbi)
  if(!all(grepl("^GO:",go_ids)) && !is.null(names(go_ids))) go_ids <- names(go_ids)
  db <- paste0('org.',species,'.eg')
  library(package=paste0(db,'.db'), character.only = T)
  eg <- mget(as.character(go_ids), get(paste0(db,'GO2ALLEGS')), ifnotfound=NA)
  x <- lapply(eg, FUN=function(x){ 
    x <- x[which(!is.na(x))]
    if(length(x)==0) return(c())
    unique(as.character(unlist(mget(as.character(x), get(paste0(db,'SYMBOL'))))))
  })
  x <- x[which(sapply(x,length)>0)]
  if(translate.ids) names(x) <- Term(names(x))
  x
}


#' findGO
#'
#' Finds GO categories where the term matches a query expression.
#'
#' @param expr The expression to query, interpreted as a regular expression if `fixed=FALSE`.
#' @param fixed Whether to interpret `expr` as a fixed expression, rather than a regular expression (default FALSE)
#' @param ontology Ontology to search, default all 3 (`BP|MF|CC`)
#'
#' @return A named vector of GO IDs, with Terms as names.
#' @export
findGO <- function(expr, fixed=FALSE, ontology="BP|MF|CC"){
  library("GO.db")
  a <- as.data.frame(GO.db::GOTERM)
  w <- intersect(grep(expr, a$Term, fixed=fixed, ignore.case=T), grep(ontology,a$Ontology))
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


#' enrichMiR2
#'
#' Looks for miRNA target enrichment across different subsets of genes
#'
#' @param genes_in_set A vector of genes in the set of interest
#' @param allgenes A vector of all genes
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param categories A list of gene lists (e.g. GO categories), with names containing the set's label
#' @param miRNA.expression A named vector of miRNAs expression values. miRNAs not in this vector are assumed to be not expressed in the system, and are not tested.
#' @param families A named vector of miRNA families, with individual miRNAs as names. If not given, internal data from the package will be used (mouse miRNA families from targetScan).
#' @param testOnlyAnnotated Whether to excluded features that are bound by no miRNA (default FALSE).
#' @param cleanNames Logical; whether to remove prefix from all miRNA names (default FALSE).
#' @param test The test to be used, either `EA`, `wEA` and `siteMir` (default `siteMir`)
#' @param minSize The minimum effective size of a category for it to be tested (default 5)
#'
#' @return a data.frame
#'
#' @export
enrichMiR2 <- function(genes_in_set, allgenes, TS, categories, miRNA.expression=NULL, families=NULL, testOnlyAnnotated=FALSE, test="siteMir", cleanNames=FALSE, minSize=5){
  test <- match.arg(test, choices=c("EA","wEA", "siteMir"))
  if(is.null(families)){
    data("miR_families")
    families <- miR_families
    if(cleanNames) names(families) <- sapply(names(families),FUN=.cleanMiRname)
  }
  if(!is.null(miRNA.expression)){
    if(is.matrix(miRNA.expression) | is.data.frame(miRNA.expression)) miRNA.expression <- rowMeans(miRNA.expression,na.rm=T)
    if(cleanNames) names(miRNA.expression) <- sapply(names(miRNA.expression),FUN=.cleanMiRname)
    miRNA.expression <- miRNA.expression[which(miRNA.expression>0)]
    families <- .filterFamilies(names(miRNA.expression), families)
    tmp <- aggregate(miRNA.expression[names(families)],by=list(family=families),na.rm=T,FUN=sum)
    fam.expr <- tmp[,2]
    names(fam.expr) <- tmp[,1]
    miRNA.expression <- list(family=fam.expr, miRNA=miRNA.expression)
  }else{
    miRNA.expression <- list(family=NULL, miRNA=NULL)
  }
  TS <- TS[which(as.character(TS$family) %in% families),]
  
  categories <- lapply(categories, y=allgenes, FUN=intersect)
  categories <- categories[which(sapply(categories,length)>=minSize)]
  res <- list()
  ll <- lapply(categories, y=genes_in_set, bg=allgenes, TS=TS, test=get(test, mode='function'), testOnlyAnnotated=testOnlyAnnotated, 
                 FUN=function(x,y,bg,TS,test,testOnlyAnnotated){
                      gin <- intersect(x,y)
                      dat <- bg %in% gin
                      names(dat) <- bg
                      m <- test(dat, TS, testOnlyAnnotated=testOnlyAnnotated)
                      m$family <- row.names(m)
                      m <- as.data.frame(m)
                      m
                 })
  catSizes <- sapply(categories, length)
  for(i in 1:length(ll)){
    ll[[i]] <- cbind(rep(names(categories)[[i]], nrow(ll[[i]])), rep(catSizes[[i]], nrow(ll[[i]])), ll[[i]])
  }
  df <- rbindlist(ll)
  colnames(df)[1:2] <- c("category","category.size")
  pc <- grep("pval",colnames(df),ignore.case=T)
  pc <- pc[length(pc)]
  df <- df[order(df$FDR,df[[pc]]),]
  return(df)
}

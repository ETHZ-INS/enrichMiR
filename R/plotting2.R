#' CDplotWrapper
#'
#' A wrapper around \code{\link{CDplot}} to work directly with DEA results and a
#' target sets as inputs.
#'
#' @param dea A named numeric vector containing the signal for each feature, or
#'   a DEA results data.frame (such as produced by major DEA packages, i.e.
#'   edgeR, DESeq2 or limma).
#' @param sets The target sets object.
#' @param setName The name of the set to plot in `sets`.
#' @param k The number of groups to form. Only applicable if by = "sites" or
#'   "score".
#' @param by The variable by which to form the groups; either "sites", "score",
#'   "best_stype" or "type". If unspecified, will determine the best depending
#'   on the available annotation data.
#' @param sameFreq Logical; whether to divide so as to obtained groups with
#'   similar frequencies (rather than groups of similar width)
#' @param pvals Logical; whether to print the p-values of KS tests between sets
#' @param line.size Size of the line.
#' @param point.size Size of the points.
#' @param ... Any other argument passed to \code{\link{CDplot}}.
#'
#' @return A ggplot.
#' @export
CDplotWrapper <- function(dea, sets, setName, k=3, 
                          by=c("auto","sites","score","best_stype","type"), 
                          sameFreq=NULL, line.size=1.2, point.size=0.8, 
                          checkSynonyms=TRUE, pvals=FALSE, ...){
  message("Preparing inputs...")
  dea <- .homogenizeDEA(dea)
  if(checkSynonyms) dea <- .applySynonyms(dea, sets)
  sets <- .list2DF(sets)
  sets <- sets[sets$set==setName,]
  if(is.null(sets$sites)) sets$sites <- 1
  by <- match.arg(by)
  if(by=="auto")
    by <- head(intersect(c("best_stype", "score", "sites"), colnames(sets)),1)
  if(by != "type" && !(by %in% colnames(sets))) 
    stop("`by` not available in `sets`.")
  if(nrow(sets)==0) stop("setName not found in `sets`.")
  if(is.null(dim(dea))){
    if(!is.numeric(dea) || is.null(names(dea)))
      stop("`dea` should be a named numeric vector or a data.frame containing ",
           "the results of a differential expression analysis.")
  }else{
    dea <- .dea2sig(.homogenizeDEA(dea), "logFC")
  }
  if(!any(sets$feature %in% names(dea)))
    stop("There seems to be no overlap between the rows of `dea` and the ",
         "features of `sets`.")
  if(is.null(sameFreq)) sameFreq <- by=="score"
  
  if(by=="best_stype" | by == "type"){
    if(by == "best_stype"){
      by <- as.character(sets[["best_stype"]])
      names(by) <- sets[["feature"]]
      by2 <- by[names(dea)]
      by2[is.na(by2)] <- "no site"
    }else{
      #get set info
      sets <- as.data.table(sets)
      sets <- melt(sets, id.vars = c("feature"),
                    measure.vars = c("Sites_8mer", "Sites_7mer_m8", 
                                     "Sites_7mer_1a"))
      sets <- sets[sets$value != 0,]
      sets$type[sets$variable == "Sites_8mer"] <- "8mer"
      sets$type[sets$variable == "Sites_7mer_m8"] <- "7mer-m8"
      sets$type[sets$variable == "Sites_7mer_1a"] <- "7mer-1a"
      #filter for those in dea and get duplicate infos
      by <- unlist(apply(sets,1, function(df){
        w <- rep(df[["type"]],df[["value"]])
        names(w) <- rep(df[["feature"]],df[["value"]])
        w}))
      #merge with dea
      by <- by[names(by) %in% names(dea)]
      by0 <- rep("no sites",length(dea[!names(dea) %in% names(by)]))
      names(by0) <- names(dea[!names(dea) %in% names(by)])
      by <- c(by,by0)
      by_df <- data.frame(feature = names(by),type = by)
      dea_df <- data.frame(feature = names(dea),logFC = dea)
      by_df <- merge(by_df,dea_df,by = "feature",all.x = TRUE)
      dea <- by_df$logFC
      names(dea) <- by_df$feature
      by2 <- by_df$type
      names(by2) <- by_df$feature
    }
    p <- CDplot(dea, by=by2, k=length(levels(as.factor(by2))), 
                sameFreq=sameFreq, size=line.size, pvals = pvals, ...)
  }else{
    by <- rowsum(sets[[by]], sets$feature)[,1]
    by2 <- by[names(dea)]
    by2[is.na(by2)] <- 0
    if(k==2){
      ll <- split(dea, by2!=0)
      names(ll) <- c("non-targets","targets")
      p <- CDplot(ll, size=line.size, pvals = pvals, ...)
    }else{
      p <- CDplot(dea, by=by2, k=k, sameFreq=sameFreq, size=line.size, 
                  pvals=pvals, ...)
    }
  }
  if(point.size>0) p <- p + geom_point(size=point.size)
  p + xlab("logFC") + ggtitle(setName)
  
}


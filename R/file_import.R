#' importTargetScan
#'
#' import a targetScan conserved targets file (conserved site and contextScore).
#'
#' @param ts.filepath Path to the targetScan conserved targets file.
#' @param featureIdentifier The type of feature identifier to use as feature name; either "symbol", "transcript", or "geneid".
#' @param recapitalizeGenes Logical; whether to recapitalize genes to the mouse standard (default FALSE)
#'
#' @return a data.frame.
#'
#' @export
importTargetScan <- function(ts.filepath, featureIdentifier="symbol", recapitalizeGenes=FALSE){
    featureIdentifier <- match.arg(featureIdentifier, c("symbol","transcript","geneID"))
    ff <- switch( tolower(featureIdentifier),
            "symbol"="Gene.Symbol",
            "transcript"="Transcript.ID",
            "geneid"="Gene.ID")
    ts2 <- read.csv(ts.filepath,header=T,stringsAsFactors=FALSE)
    tsg <- aggregate(ts2[,c("Total.num.conserved.sites", "Cumulative.weighted.context...score")],by=list(family=ts2$miRNA.family, representative=ts2$Representative.miRNA, gene=recapitalizeGenes(ts2[[ff]])),FUN=function(x){ return(max(abs(x))) })
    colnames(tsg) <- c("family","rep.miRNA","feature","sites","score")
    return(tsg)
}

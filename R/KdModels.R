#' @exportClass KdModel
setClass(
  "KdModel",
  contains="lm",
  validity=function(object){
    if(is.null(object$name) || !is.character(object$name) ||
       length(object$name)!=1)
      stop("The model should have a `name` slot.")
    if(is.null(object$canonical.seed) || !is.character(object$canonical.seed) ||
       length(object$canonical.seed)!=1)
      stop("The model should have a `canonical.seed` slot (character).")
  }
)


#' @export
setMethod("show", "KdModel", function(object){
  cat(paste0("A `KdModel` for ", object$name, " (", object$canonical.seed,")"))
})

#' @export
setMethod("summary", "KdModel", function(object){
  co <- coefficients(object)
  co <- co[grep("^sr", names(co))]
  co <- co[grep(":",names(co),invert=TRUE)]
  names(co) <- gsub("^sr","",names(co))
  sort(co)
})

#' @export
setMethod("predict", signature(object="KdModel"), function(object, kmers){
  predictKD(kmers, object)
})

#' getKdModel
#' 
#' Summarizes the binding affinity of 12-mers using linear models.
#'
#' @param kd A data.frame with at least the columns "X12mer" and "log_kd", or the path to
#' such a data.frame
#' @param name The optional name of the miRNA
#'
#' @return A linear model of class `KdModel`
#' @export
getKdModel <- function(kd, name=NULL){
  if(is.character(kd) && length(kd)==1){
    if(is.null(name)) name <- gsub("\\.rds$","",basename(kd),ignore.case=TRUE)
    kd <- read.delim(kd, header=TRUE)
  }
  if("mirseq" %in% colnames(kd)){
    mirseq <- as.character(kd$mirseq[1])
    seed <- as.character(reverseComplement(DNAString(substr(mirseq, 2,8))))
    name <- as.character(kd$mir[1])
    kd <- kd[,c("X12mer","log_kd")]
  }else{
    mirseq <- seed <- NULL
  }
  kd <- kd[grep("X",kd$X12mer,invert=TRUE),]
  pwm <- Biostrings::consensusMatrix(
    as.character(rep(kd$X12mer, floor( (10^(-kd$log_kd))/10 ))),
    as.prob=TRUE
  )
  fields <- c("sr","A","fl")
  if(!all(fields %in% colnames(kd)))
    kd <- prep12mers(kd)
  fields <- c(fields, "log_kd")
  if(!all(fields %in% colnames(kd))) stop("Malformed `kd` data.frame.")
  w <- (1-kd$log_kd)^2
  w[which(w<0.5)] <- 0.5
  mod <- lm( log_kd~sr*A+fl, data=kd, model=FALSE, x=FALSE, y=FALSE )
  mod$cor.with.cnn <- cor(mod$fitted.values, kd$log_kd)
  mod$mae.with.cnn <- median(abs(mod$residuals))
  mod$residuals <- NULL
  mod$fitted.values <- NULL
  mod$weights <- NULL
  mod$assign <- NULL
  mod$effects <- NULL
  mod$qr <- list(pivot=mod$qr$pivot)
  mod$name <- name
  mod$mirseq <- mirseq
  mod$canonical.seed <- seed
  mod$pwm <- pwm
  class(mod) <- c("KdModel", class(mod))
  new("KdModel", mod)
}

#' prep12mers
#'
#' @param x A vector of 12-mers, or a data.frame containing at least the columns
#' "X12mer" and "log_kd"
#' @param mod An optional linear model summarizing the kd activity.
#' @param maxSeedMedian Max median log_kd for alternative seed inclusion.
#' @param maxNSeeds Maximum number of seeds to include in the model
#'
#' @return A data.frame
#' @export
prep12mers <- function(x, mod=NULL, maxSeedMedian=-1.2, maxNSeeds=30){
  if(is.data.frame(x)){
    if(!all(c("X12mer","log_kd") %in% colnames(x)))
      stop("`x` should be a character vector or a data.frame with the columns ",
           "'X12mer' and 'log_kd'")
    x <- x[grep("X",x$X12mer,invert=TRUE),]
    x <- cbind( x[,intersect(c("log_kd","energy"), colnames(x)),drop=FALSE], 
                prep12mers(x$X12mer) )
    ag <- aggregate(x$log_kd, by=list(seed=x$sr), FUN=median)
    seedMed <- ag$x
    names(seedMed) <- ag[,1]
    seedMed <- sort(seedMed[seedMed<=maxSeedMedian])
    seedMed <- names(seedMed)[seq_len(min(length(seedMed), maxNSeeds))]
    d <- data.frame( log_kd=0, sr="other", A=FALSE, fl=levels(x$fl)[1] )
    if("energy" %in% colnames(x)) d$energy <- 0
    return( rbind( d[,colnames(x)], x[x$sr %in% seedMed,] ) )
  }
  x <- as.character(x)
  sr <- sapply(x, FUN=function(x) substr(x, 3,9))
  if(!is.null(mod)){
    if(is.null(mod$xlevels$sr)) stop("The model contains no seed levels.")
    sr.lvls <- mod$xlevels$sr
    sr[!(sr %in% sr.lvls)] <- "other"
  }else{
    sr.lvls <- c("other",unique(sr))
  }
  sr <- factor(sr, sr.lvls)
  fl <- paste0(substr(x,1,2), substr(x,11,12))
  d <- data.frame( sr=sr, A=as.logical(substr(x, 10, 10)=="A"),
                   fl=factor(fl, levels=getKmers(4)), row.names=NULL )
  d$A[is.na(d$A)] <- FALSE
  if(!is.null(mod) && any(is.na(d$fl))){
    fl <- sort(coef(mod)[grep("fl",names(coef(mod)))])
    d$fl[is.na(d$fl)] <- gsub("^fl","",names(fl)[floor(length(fl)/2)])
  }
  d
}

#' getKmers
#'
#' Returns all combinations of `n` elements of `from`
#'
#' @param n Number of elements
#' @param from Letters sampled
#'
#' @return A character vector
#' @export
#'
#' @examples
#' getKmers(3)
getKmers <- function(n=4, from=c("G", "C", "T", "A")){
  apply(expand.grid(lapply(seq_len(n), FUN=function(x) from)),
        1,collapse="",FUN=paste)
}


#' plotKdModel
#'
#' @param mod A `KdModel`
#' @param what Either 'seeds', 'logo', or 'both' (default)
#'
#' @return A plot
#' @export
plotKdModel <- function(mod, what=c("both","seeds","logo")){
  library(ggplot2)
  what <- match.arg(what)
  if(what=="seeds"){
    coe <- coefficients(mod)
    medfl <- median(coe[grep("^fl",names(coe),value=TRUE)],na.rm=TRUE)
    co <- -coe["(Intercept)"]-coe[paste0("sr",mod$xlevels$sr[-1])]-medfl
    names(co) <- gsub("^sr","",names(co))
    co <- sort(co)
    co <- data.frame(seed=factor(names(co), names(co)), log_kd=as.numeric(co))
    co$type <- sapply(as.character(co$seed), seed=mod$canonical.seed, .getMatchType)
    coA <- co
    aint <- coe[paste0("sr",coA$seed,":ATRUE")]
    aint[is.na(aint)] <- 0
    coA$log_kd <- - coe["ATRUE"] - aint
    coA$type <- "+A"
    co <- rbind(co,coA)
    co$type <- factor(co$type, c("+A","7mer-m8","6mer","offset 6mer","non-canonical"))
    p <- ggplot(co, aes(seed, log_kd, fill=type)) + geom_col() + 
      coord_flip() + ylab("-log_kd")
    if(!is.null(mod$name)) p <- p + ggtitle(mod$name)
    return( p )
  }
  if(what=="logo") return(seqLogo::seqLogo(mod$pwm, xfontsize=12, yfontsize=12))
  cowplot::plot_grid( plotKdModel(mod, "seeds"),
                       grid::grid.grabExpr(plotKdModel(mod, "logo")),
                       nrow=2, rel_heights=c(6,4))
}

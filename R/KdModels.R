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
  mod <- lm( log_kd~sr*A+fl, data=kd, model=FALSE, weights=w, x=FALSE, y=FALSE )
  mod$cor.with.cnn <- cor(mod$fitted.values, kd$log_kd)
  mod$mae.with.cnn <- median(abs(mod$fitted.values-kd$log_kd))
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
  mod
}

#' prep12mers
#'
#' @param x A vector of 12-mers, or a data.frame containing at least the columns
#' "X12mer" and "log_kd"
#' @param mod An optional linear model summarizing the kd activity.
#' @param maxSeedMedian Max median log_kd for alternative seed inclusion.
#'
#' @return A data.frame
#' @export
prep12mers <- function(x, mod=NULL, maxSeedMedian=-1.2){
  if(is.data.frame(x)){
    if(!all(c("X12mer","log_kd") %in% colnames(x)))
      stop("`x` should be a character vector or a data.frame with the columns ",
           "'X12mer' and 'log_kd'")
    x <- x[grep("X",x$X12mer,invert=TRUE),]
    x <- cbind(x[,"log_kd",drop=FALSE], prep12mers(x$X12mer))
    ag <- aggregate(x$log_kd, by=list(seed=x$sr), FUN=median)
    ag <- as.character(ag[ag$x <= maxSeedMedian,1])
    return( rbind( data.frame( log_kd=0, sr="other", A=FALSE, 
                               fl=levels(x$fl)[1] ),
                   x[x$sr %in% ag,] ) )
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
  fl <- sort(coef(mod)[grep("fl",names(coef(mod)))])
  d$fl[is.na(d$fl)] <- gsub("^fl","",names(fl)[floor(length(fl)/2)])
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
    # mer8 <- co[nrow(co),,drop=FALSE]
    # mer8$log_kd <- mer8$log_kd-coe[paste0("sr",mer8$seed,":ATRUE")]-coe["ATRUE"]
    # mer8$type="8mer"
    # co <- rbind(co, mer8)
    # co$type <- factor(co$type, c("8mer","7mer-m8","6mer","non-canonical"))
    p <- ggplot(co, aes(seed, log_kd, fill=type)) + geom_col() + 
      coord_flip() + ylab("-log_kd")
    if(!is.null(mod$name)) p <- p + ggtitle(mod$name)
    return( p )
  }
  if(what=="logo") return(seqLogo::seqLogo(mod$pwm))
  library(cowplot)
  plot_grid( plotKdModel(mod, "seeds"),
             grid::grid.grabExpr(plotKdModel(mod, "logo")),
             nrow=2)
}


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
    #sr <- c("other",gsub("^sr","",names(co)[seq.int(2,col)]))
    data.frame( seed=names(co)[seq_len(col)],
                seed.coef=as.integer(round(100*co[seq_len(col)])),
                seed.A=as.integer(round(100*co[seq.int(col+1,length(co))])),
                stringsAsFactors = FALSE
                )
  }), .id="miRNA")
  sr$miRNA <- factor(sr$miRNA)

  otherfields <- c( "rank","qr","df.residual","mirseq","canonical.seed","pwm",
                    "cor.with.cnn","mae.with.cnn","name" )
  names(otherfields) <- otherfields
  other <- lapply(mods, FUN=function(x){
    y <- lapply(otherfields, FUN=function(f) x[[f]])
    y$srlvls <- x$xlevels$sr
    y
  })
  
  modf <- mods[[1]]
  modf$xlevels <- list(sr=c(), fl=modf$xlevels$fl)
  for(f in otherfields) modf[[f]] <- NULL
  
  mods <- list(frame=modf, sr=sr, fl=fl, other=other)
  class(mods) <- c("CompressedKdModelList", "list")
  mods
}


decompressKdModList <- function(mods){
  if(length(mods)==1 && is.character(mods)) mods <- readRDS(mods)
  if(!is(mods,"CompressedKdModelList")) stop("`mods` is not a CompressedKdModelList")
  otherfields <- c( "rank","qr","df.residual","mirseq","canonical.seed","pwm",
                    "cor.with.cnn","mae.with.cnn","name" )
  SR <- split(mods$sr, mods$sr$miRNA)
  names(nn) <- nn <- names(mods$other)
  lapply(nn, FUN=function(n){
    mod <- mods$frame
    for(f in otherfields){
      if(f %in% names(mods$other[[n]])) mod[[f]] <- mods$other[[n]][[f]]
    }
    mod$xlevels$sr <- mods$other[[n]]$srlvls
    co <- SR[[n]]
    co2 <- c(co[,3],co[,4])/100
    names(co2) <- c(co$seed,paste0(co$seed,":ATRUE"))
    names(co2)[nrow(co)+1] <- "ATRUE"
    fl <- mods$fl[,n]/100
    names(fl) <- row.names(mods$fl)
    mod$coefficients <- c( co2[grep(":",names(co2),invert=TRUE)],
                           fl, co2[grep(":",names(co2))] )
    mod
  })
}

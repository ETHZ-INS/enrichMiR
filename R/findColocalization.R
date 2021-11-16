#' testColocalization
#'
#' Creates an enrich.results object and performs enrichment analysis of miRNA
#' colocalizations.
#'
#' In detail, finds colocalizations of miRNAs with other miRNAs or user-defined
#' Partners that are located on the same gene or transcript within a specified
#' distance. Fisher's enrichment test is applied in order to determine which
#' miRNA-pairs are significantly enriched in the particular universe of
#' expressed genes.
#'
#' In "global" mode, statistical enrichment analysis is performed on the
#' universe of expressed genes, meaning is there a statistical enrichment of a
#' certain "pair" with regards to the other pairs that are found in the
#' universe. The results of this search can be accessed with "Fisher.Pair" in
#' the enrich.result object.
#'
#' In "subset" mode, we search for enrichments of pairs in either specified
#' genes of interest or up- / downregulated genes of a differential expression
#' analysis. By default, three independent statistics will be performed: 1)
#' "Fisher.Pair": In order to check whether there is an enrichment of a specific
#' pair in the subset vs the background. 2) "Fisher.MIR.Co": Check whether the
#' specified microRNA has more colocolizations with this partner than with any
#' other. 3) "Fisher.Partner.Co": Check whether the specified Partner has more
#' colocolizations with this microRNA than with any other.
#'
#'
#'
#' @param x The signature in which to look for a signal. This can be either:
#'   \itemize{ \item A data.frame of the results of a differential expression
#'   analysis, with features (e.g. genes) as row names and with at least the
#'   following columns: `logFC`, `FDR`; \item A logical vector of membership to
#'   the geneset of interest, with feature (e.g. genes) as names, in which case
#'   only tests based on a binary signal will be available. \item A vector of
#'   feature names; in this case the search mode is "global" and only tests
#'   based on a binary signal will be available. } By convention, we only test
#'   for enrichments in annotated genes.
#' @param background A character vector of background; ignored if `x` is not a
#'   character vector.
#' @param mir_pos A data.frame (or DataFrame) with positional information of
#'   miRNA Binding Sites as obtained through the findSeedMatches function or
#'   enrichMiR:::.getTargetScanPos. (Here one should the correct link) Should
#'   contain at least the following columns: `set`, `feature`,`start` & `stop`.
#' @param sets_pos A data.frame (or DataFrame) / data.table with positional
#'   information of any RBP motif sites as for example obtained through the
#'   ""link to RBP_Pos function"" or in custom manner. Should contain at least
#'   the following columns: `set` = RBP or Name, `feature` = Gene or Transcript
#'   Target,`start` & `stop` (positional information, with regards to the
#'   3'UTR).
#' @param mir_pos.properties Any further information about the mir_pos object;
#'   this can either be a data.frame (or DataFrame), with row.names
#'   corresponding to names of `sets` (or to alternative names), or a named
#'   vector (e.g. miRNA expression values). If this is given, sets (in the
#'   mir_pos object) not included within this object will be discarded.
#' @param th.abs.logFC The minimum absolute log2-foldchange threshold for a
#'   feature to be considered differentially-expressed (default 0). Ignored if
#'   `x` is not a DEA data.frame.
#' @param th.FDR The maximum FDR for a feature to be considered differentially-
#'   expressed (default 0.05). Ignored if `x` is not a DEA data.frame.
#' @param minSize The minimum size of a set to be tested (default 5).
#' @param maxSize The maximum size of a set to be tested (default Inf).
#' @param co_dist_min Minimum number of nucleotides to be located between two
#'   miRNA-seeds or a miRNA-seed and RBP-motif. By specifying co_dist_min and
#'   co_dist_max, the user can search for certain types of functional miRNA
#'   colocalization. Mutual miRNA blocking is suggested to happen within a
#'   distance of 0-7Nt. Cooperative repression of two miRNAs is reported to
#'   happen within a distance of 8-40Nt, and at bigger distance, miRNAs are
#'   suggested to act independently on given targets. If co_dist_max = NULL,
#'   enrichment analysis will be performed on the whole transcript.
#' @param co_dist_max Maximum number of nucleotides to be located between two
#'   miRNA-seeds or a miRNA-seed and RBP-motif. For further info, see
#'   co_dist_min.
#' @param force.global.mode To force a statistical enrichment analysis on the
#'   whole universe of expressed genes, even in case there is for example a
#'   data.frame of the results of a differential expression analysis, specified
#'   as x. By default: FALSE.
#' @param BPPARAM could be potentially included, see .findcolocPL at the bottom
#'   \link{BiocParallel} multithreading parameters. Used to multithread
#'   colocalization.
#'   
#' @return an enrich.results object.
#'
#' @import S4Vectors
#' @importFrom BiocParallel bplapply SerialParam
#' @export
testColocalization <- function(x, background=NULL, mir_pos, sets_pos = NULL,
                            mir_pos.properties=NULL, th.abs.logFC=0, th.FDR=0.05, 
                            minSize=5, maxSize=Inf, co_dist_min = 0, co_dist_max = NULL, 
                            force.global.mode = FALSE, BP=NULL, 
                            ...){
  
  
  if(is.null(co_dist_min) && is.null(co_dist_max))
    stop("Either a minium or maximum colocalization distance has to be chosen. If you want to get all pairs, set 'co_dist_min = 0'")
  
  if(!is.null(co_dist_max)){
  if(co_dist_min != 0 && co_dist_min >= co_dist_max)
    stop("The specified maximal distance between two motifs has to be bigger than the minimal distance")
  }
  
  ## Data Preparation
  
  mode <- "subset"

  if(is.character(x)){
    if(is.null(background)) {
        message("Search for colocalization of motifs within the full given set of genes")
        mode <- "global"
        t <- !x %in% x
        names(t) <- x
        x <- t
      }else{
        if(force.global.mode){
          mode = "global"
          t <- rep(FALSE,length(x))
          names(t) <- x
          x <- t
          warning("`background` ignored.")
        }else{
          x <- background %in% x
          names(x) <- background 
        }
        
      }
  }
  
  
  if(!is.null(dim(x))){ 
    x <- .homogenizeDEA(x)
    if(force.global.mode){
      mode = "global"
      t <- !row.names(x) %in% row.names(x)
      names(t) <- row.names(x)
      x <- t
    }}
  

  x <- .applySynonyms(x, mir_pos)
  
  if( (is.null(dim(x)) && (!any(names(x) %in% mir_pos$feature))) ||
      (!is.null(dim(x)) && (!any(row.names(x) %in% mir_pos$feature))) )
    stop("There is no match between the features of the `mir_pos-object` and those of `x`. ",
         "Are you sure that you are using the right annotation for your data?")
    
  ## Prepare the MIRS-Object
  
  # restrict the mir_pos object according to properties and sizes
  mir_pos.properties <- .filterMatchSets(mir_pos, mir_pos.properties)
  mir_pos <- mir_pos[mir_pos$set %in% row.names(mir_pos.properties),]
  mir_pos.properties$set_size <- table(mir_pos$set)[row.names(mir_pos.properties)]
  mir_pos.properties <- mir_pos.properties[mir_pos.properties$set_size >= minSize & 
                                             mir_pos.properties$set_size <= maxSize,
                                     , drop=FALSE]
  mir_pos <- mir_pos[mir_pos$set %in% row.names(mir_pos.properties),]
  if(is.factor(mir_pos$set)) mir_pos$set <- droplevels(mir_pos$set)
  
  
  # filter the mir_pos object for expressed genes
  if(!is.null(x)){
    if(is.null(dim(x))){
      mir_pos <- mir_pos[mir_pos$feature %in% names(x),]
    }else{
      mir_pos <- mir_pos[mir_pos$feature %in% row.names(x),]
    }
  }
  if(length(mir_pos$feature) == 0)
    stop("No miRNA pairs found in expressed genes")
  
  
  # Create the signatures
  signature.list <- list()
  if(mode == "subset" && !is.null(dim(x))){
    signature.list <- list(
      down=.dea2binary(x, th=th.FDR, th.alfc=th.abs.logFC, restrictSign=-1, ...),
      up=.dea2binary(x, th=th.FDR, th.alfc=th.abs.logFC, restrictSign=1, ...)
    )
  }else{
    signature.list <- list(sig = x)
  }
  
  ll <- lapply(signature.list,FUN=function(x) names(x)[x])
  if(length(ll) > 1) {
    mir_pos$signal_info <- ifelse(mir_pos$feature %in% ll$up,"up",ifelse(mir_pos$feature %in% ll$down,"down","background"))
  }else{
    mir_pos$signal_info <- ifelse(mir_pos$feature %in% ll$sig,"sig","background")
  }
  
  #Prepare the GRanges Object
  mirs <- GRanges(seqnames=mir_pos$feature, IRanges(mir_pos$start, mir_pos$end), miRNA = mir_pos$set, signal = mir_pos$signal_info)
  mirs.back <- mirs[mirs$signal == "background",]
  if(mode=="subset"){
    if(!is.null(dim(x))){
      mirs.up <- mirs[mirs$signal == "up",]
      mirs.down <- mirs[mirs$signal == "down",]
    }else{
      mirs.sig <- mirs[mirs$signal == "sig",]
    }
  }
  
  ## Potentially create the RBP-Object
  
  # filtering etc.
  if(!is.null(sets_pos)){
    if(is.null(dim(sets_pos)))
      stop("sets_pos is not in the correct format, see description")
    sets_pos <- DataFrame(sets_pos)
    
    # Get the correct names
    ff <- split(sets_pos$feature,sets_pos$feature)
    ff <- data.frame(row.names=names(ff), feature=names(ff),stringsAsFactors = FALSE)
    
    if(!is.null(metadata(mir_pos)$feature.synonyms)) {
      ss <- metadata(mir_pos)$feature.synonyms
      if(any(row.names(ff) %in% names(ss))){
        ff <- ff[order(ff$feature),,drop = FALSE]
        n <- row.names(ff)
        w <- which(n %in% names(ss))
        n[w] <- as.character(ss[n[w]])
        w <- duplicated(n)
        ff <- ff[!w,,drop = FALSE]
        row.names(ff) <- n[!w]
      }}
    w <- which(row.names(ff) %in% mir_pos$feature)
    ff <- ff[w,,drop = FALSE]  
    sets_pos <- sets_pos[sets_pos$feature %in% ff$feature,]
    colnames(ff)[which(colnames(ff) == "feature")] <- "merge_name"
    ff$feature <- row.names(ff)
    colnames(sets_pos)[which(colnames(sets_pos) == "feature")] <- "merge_name"
    sets_pos <- merge(sets_pos,ff,by = "merge_name",all.x = TRUE)
    
    #filter for the set.size
    sets_pos.properties <- .filterMatchSets(sets_pos,props = NULL)
    sets_pos.properties$set_size <- table(sets_pos$set)[row.names(sets_pos.properties)]
    sets_pos.properties <- sets_pos.properties[sets_pos.properties$set_size >= minSize & 
                                                sets_pos.properties$set_size <= maxSize,
                                              , drop=FALSE]
    sets_pos <- sets_pos[sets_pos$set %in% row.names(sets_pos.properties),]
    if(is.factor(sets_pos$set)) sets_pos$set <- droplevels(sets_pos$set)
    
    # Get the signal info
    if(mode == "subset" && !is.null(dim(x))){
      sets_pos$signal_info <- ifelse(sets_pos$feature %in% ll$up,"up",ifelse(sets_pos$feature %in% ll$down,"down","background"))
    }else{
      sets_pos$signal_info <- ifelse(sets_pos$feature %in% ll$sig,"sig","background")
    }
    
    
    #Prepare the GRanges Object
    rbps <- GRanges(seqnames=sets_pos$feature, IRanges(sets_pos$start, sets_pos$end), RBP = sets_pos$set, signal = sets_pos$signal_info)
    rbps.back <- rbps[rbps$signal == "background",]
    if(mode=="subset"){
      if(!is.null(dim(x))){
        rbps.up <- rbps[rbps$signal == "up",]
        rbps.down <- rbps[rbps$signal == "down",]
      }else{
        rbps.sig <- rbps[rbps$signal == "sig",]
      }
    }
    
  }
  ## I guess here one could insert some memory optimization steps
   
  
  ## Create the Colocalization Tables
  
  if(is.null(BP)) BP <- SerialParam()
  # global mode
  # all genes should be assigned to background in global mode
  if(mode == "global"){
    if(!is.null(sets_pos)){
      dt <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs.back,rbps.back,BP) 
    }else{
      dt <- .findcoloc1mir(co_dist_min,co_dist_max,mirs.back,BP)
    }
  }else{
    # subset mode
    # prepare the background table
    if(!is.null(sets_pos)){
      dt.back <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs.back, rbps.back,BP)
    }else{
      dt.back <- .findcoloc1mir(co_dist_min,co_dist_max,mirs.back,BP)
    }
    colnames(dt.back)[-c(1,2)] <- paste(colnames(dt.back)[-c(1,2)], "in.back", sep = "_")
        
    # prepare the subset tables and merge
    if(!is.null(dim(x))){
      if(!is.null(sets_pos)){
        dt.up <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs.up,rbps.up,BP)
        dt.down <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs.down, rbps.down,BP)
      }else{
        dt.up <- .findcoloc1mir(co_dist_min,co_dist_max,mirs.up,BP)
        dt.down <- .findcoloc1mir(co_dist_min,co_dist_max,mirs.down,BP)
      }
      colnames(dt.up)[-c(1,2)] <- paste(colnames(dt.up)[-c(1,2)], "in.up", sep = "_")
      colnames(dt.down)[-c(1,2)] <- paste(colnames(dt.down)[-c(1,2)], "in.down", sep = "_")
      dt <- merge(dt.up,dt.back,by = c("miRNA","Partner"), all = TRUE)
      dt <- merge(dt.down,dt,by = c("miRNA","Partner"), all = TRUE)
      }else{
        if(!is.null(sets_pos)){
           dt.sig <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs.sig, rbps.sig,BP)
        }else{
          dt.sig <- .findcoloc1mir(co_dist_min,co_dist_max,mirs.sig,BP)
        }
        colnames(dt.sig)[-c(1,2)] <- paste(colnames(dt.sig)[-c(1,2)], "in.sig", sep = "_")
        dt <- merge(dt.sig,dt.back,by = c("miRNA","Partner"), all = TRUE)
      }
    dt[is.na(dt)] <- 0
    }
    
  ## Create the Enrich.Result Object
  
  
  if(!is.null(sets_pos)){
    o <- new("enrich.results", input=list(x=x, sets.properties= 
                                          list("miRNA" = mir_pos.properties, "Partner" = sets_pos.properties)), 
           binary.signatures=signature.list, info=list(call=match.call(),type = "colocalization"))
  }else{
    o <- new("enrich.results", input=list(x=x, sets.properties= 
                                            list("miRNA" = mir_pos.properties)), 
             binary.signatures=signature.list, info=list(call=match.call(),type = "colocalization"))
  }
  
  # add mir-overlap info
  if(mode == "subset"){
    o@overlaps <- lapply(signature.list, FUN=function(x){
    FactorList(lapply(split(mir_pos$feature, mir_pos$set), 
                      FUN=function(y) sort(intersect(y,names(x)[x]))))
  })
  }
  if(mode == "global"){
    o@overlaps <- list("overlaps" = FactorList(lapply(split(mir_pos$feature, mir_pos$set), 
                        FUN=function(y) sort(intersect(y,names(x)))))
    )
  }
  
  
  ## Now comes the statistic
  
  # global mode
  if(mode == "global"){
    dt <- .fisher_test(dt,mod = "global")
    o@res["Fisher.Pair"] <- list(dt)
  }
  

  # subset mode
  if(mode == "subset"){
    if(!is.null(dim(x))){
      
      # 1) Check for the enrichment of a specific pair in the subset vs the background.
      dt.up <- dt[,c(1,2,7,10,11,14)]
      dt.up <- .fisher_test(dt.up,mod = "subset")
      o@res["Fisher.Pair.Up"] <-list(dt.up)
      dt.down <- dt[,c(1,2,3,6,11,14)]
      dt.down <- .fisher_test(dt.down,mod = "subset")
      o@res["Fisher.Pair.Down"] <-list(dt.down)
    
      # 2) Check whether the specified microRNA has more colocolizations with this partner than
      #    expected by chance.
      dt.up <- dt[,c(1,2,7,8,11,12)]
      dt.up <- .fisher_test(dt.up,mod = "subset")
      o@res["Fisher.MIR.Co.Up"] <-list(dt.up)
      dt.down <- dt[,c(1,2,3,4,11,12)]
      dt.down <- .fisher_test(dt.down,mod = "subset")
      o@res["Fisher.MIR.Co.Down"] <-list(dt.down)
      
      
      if(!is.null(sets_pos)){
      # 3) Check whether the specified Partner has more colocolizations with this microRNA than
      #    expected by chance.
      dt.up <- dt[,c(1,2,7,9,11,13)]
      dt.up <- .fisher_test(dt.up,mod = "subset")
      o@res["Fisher.Partner.Co.Up"] <-list(dt.up)
      dt.down <- dt[,c(1,2,3,5,11,13)]
      dt.down <- .fisher_test(dt.down,mod = "subset")
      o@res["Fisher.Partner.Co.Down"] <-list(dt.down)
      }
      
    }else{
      
      # 1) Check for the enrichment of a specific pair in the subset vs the background.
      dt.sig <- dt[,c(1,2,3,6,7,10)]
      dt.sig <- .fisher_test(dt.sig,mod = "subset")
      o@res["Fisher.Pair"] <-list(dt.sig)
      
      # 2) Check whether the specified microRNA has more colocolizations with this partner than
      #    expected by chance.
      dt.sig <- dt[,c(1,2,3,4,7,8)]
      dt.sig <- .fisher_test(dt.sig,mod = "subset")
      o@res["Fisher.MIR.Co"] <-list(dt.sig)
      
      if(!is.null(sets_pos)){
      # 3) Check whether the specified Partner has more colocolizations with this microRNA than
      #    expected by chance.
      dt.sig <- dt[,c(1,2,3,5,7,9)]
      dt.sig <- .fisher_test(dt.sig,mod = "subset")
      o@res["Fisher.Partner.Co"] <-list(dt.sig)
      }
      
    } 
  }
  
  o
}

  

## Helper Functions

  
.fisher_test <- function(dt,mod){
  dt$pvalue <- apply(dt,1,function(dt){
    m <- matrix(as.numeric(dt[c(3, 4, 5, 6)]), ncol = 2)
    f <- fisher.test(as.table(m), alt="greater")
    return(f$p.value)
  })
  dt$FDR <- p.adjust(dt$pvalue, method = "BH")
  dt <- as.data.frame(dt)
  if(mod == "subset"){
    dt$enrichment <- apply(dt,1,function(dt){
      c <- as.numeric(dt[c(3, 4, 5, 6)])
      c <- ((c[1] + 0.25)/max(c[2],1))/((c[3]+0.25)/max(c[4],1))
      ifelse(c==0,0,log2(c))
    })
  }
  dt
}  




.findcoloc1mir <- function(min.dist,max.dist,Pos_object,BP) {
  library(BiocParallel)
  if(is.null(BP)) BP <- SerialParam()
  
  df <- sapply( split(Pos_object, Pos_object$miRNA), FUN=function(mirs_all){
    # each mirs contains all the binding sites of a single miRNA
    # this means basically that a second argument (the mirs_all) is passed to the lapply function
    ll <- bplapply(split(Pos_object, Pos_object$miRNA), mi=mirs_all,BPPARAM=BP, FUN=function(x,mi){
      suppressWarnings(distanceToNearest(x, mi, ignore.strand=T)@elementMetadata$distance)
    })
    sapply(ll,FUN=function(x) {
      #assign the min-max values
      if(is.numeric(min.dist)){
        if(is.numeric(max.dist)){
          sum(abs(x)>= min.dist & abs(x)<= max.dist, na.rm = T)
        }else{
          sum(abs(x)>= min.dist, na.rm = T)
        }
      }else{
        if(is.numeric(max.dist)){
          sum(abs(x)<= max.dist, na.rm = T)  
        }else{
          stop("This stop actually shouldn't happen")
        }
      }
    })
  })
  #Prepare a dataframe to analyze the results and filter out the "self-counts" of each micro which are on the diagonale of the df.
  df<- as.data.frame(df)
  df[upper.tri(df, diag = TRUE)] <- NA
  df$miRNA <- row.names(df)
  dt <- as.data.table(df)
  dt <- melt(dt, id.vars = 'miRNA', na.rm = TRUE, variable.name = "Partner", value.name = "overlap")
  
  # sum up the pairs
  mir_nr <- dt[,.(sum_mir = sum(.SD)),by = "miRNA", .SDcols = c("overlap")]
  partner_nr <- dt[,.(sum_partner = sum(.SD)),by = "Partner", .SDcols = c("overlap")]
  pairs <- merge(mir_nr,partner_nr, by.x = 'miRNA', by.y = 'Partner', all=TRUE)
  pairs[is.na(pairs)] <- 0
  pairs$nr_of_pairs <- pairs$sum_mir + pairs$sum_partner
  
  # get miRNA numbers for statistics
  dt <- merge(dt,pairs[,c("miRNA","nr_of_pairs")], by = 'miRNA', all.x = TRUE)
  colnames(dt)[colnames(dt)=="nr_of_pairs"] <- "nr_of_pairs_miRNA"
  dt <- merge(dt,pairs[,c("miRNA","nr_of_pairs")], by.x = 'Partner', by.y = 'miRNA', all.x = TRUE)
  colnames(dt)[colnames(dt)=="nr_of_pairs"] <- "nr_of_pairs_Partner"
  dt$other_pairs_miRNA <- dt$nr_of_pairs_miRNA - dt$overlap
  dt$other_pairs_Partner <- dt$nr_of_pairs_Partner - dt$overlap
  dt$background_pairs <- sum(dt$overlap) - (dt$overlap + dt$other_pairs_miRNA + dt$other_pairs_Partner)
  dt <- dt[,-c(4,5)]
  dt <- setcolorder(dt,neworder = c("miRNA","Partner"))
  dt
}





.findcoloc1rbp<- function(min.dist,max.dist,mirs,rbps,BP) {
  library(BiocParallel)
  if(is.null(BP)) BP <- SerialParam()
  # This should be miRNA centered, correct? Would that matter?
  
  df <- sapply( split(mirs, mirs$miRNA), FUN=function(mirs_all){
    # mirs contains all the binding sites of a single miRNA
    ll <- bplapply(split(rbps, rbps$RBP), mi=mirs_all, BPPARAM=BP,FUN=function(x,mi){
      suppressWarnings(distanceToNearest(x, mi, ignore.strand=T)@elementMetadata$distance)
    })
    sapply(ll,FUN=function(x) {
      #assign the min-max values
      if(is.numeric(min.dist)){
        if(is.numeric(max.dist)){
          sum(abs(x)>= min.dist & abs(x)<= max.dist, na.rm = T)
        }else{
          sum(abs(x)>= min.dist, na.rm = T)
        }
      }else{
        if(is.numeric(max.dist)){
          sum(abs(x)<= max.dist, na.rm = T)  
        }else{
          stop("This stop actually shouldn't happen")
        }
      }
    })
  })
  
  dt <- as.data.table(df, keep.rownames = TRUE)
  colnames(dt)[colnames(dt)=="rn"] <- "Partner"
  dt <- melt(dt, id.vars = 'Partner', na.rm = TRUE, variable.name = "miRNA", value.name = "overlap")
  
  # sum up the pairs
  mir_nr <- dt[,.(sum_mir = sum(.SD)),by = "miRNA", .SDcols = c("overlap")]
  partner_nr <- dt[,.(sum_partner = sum(.SD)),by = "Partner", .SDcols = c("overlap")]
  
  # get miRNA numbers for statistics
  dt <- merge(dt,mir_nr, by = 'miRNA', all.x = TRUE)
  dt <- merge(dt,partner_nr, by = 'Partner', all.x = TRUE)
  dt$other_pairs_miRNA <- dt$sum_mir - dt$overlap
  dt$other_pairs_Partner <- dt$sum_partner - dt$overlap
  dt$background_pairs <- sum(dt$overlap) - (dt$overlap + dt$other_pairs_miRNA + dt$other_pairs_Partner)
  dt <- dt[,-c(4,5)]
  dt <- setcolorder(dt,neworder = c("miRNA","Partner"))
  dt
  
}





  
  
  
  
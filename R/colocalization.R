

# split_by seed in TargetScan-mode
# split_by miRNA in KD-mode

# co_dist between seeds in TS
# co_dist between 12mers in KD

# For the KDs, include option to work on a Transcript level?

#0) ?berlege Dir 1-2 leichte Analysen um sie Gerhard und den anderen zu zeigen (Syncrip and PTX RBPs + Celf from Postar?)
#   Human Synaptogenesis genes, enrichment for RBP + miRNA

#3) test ob das mit dem custom funktioniert (>> vll direkt mit Syncrip um es Gerhard zu zeigen?)

#5) Prepare all the RBPPos objects for each species 






#' Find miRNA colocalizations
#'
#' Finds pairs of miRNAs that are located on the same gene or transcript within
#' a distance that the user can specify. Fisher's enrichment test is applied in
#' order to determine which miRNA-pairs are significantly enriched in the
#' particular universe of expressed genes.
#'
#' @param expressed.genes A vector of expressed genes (at the moment EnsemblID,
#'   maybe give option?) as a universe to search in
#' @param search.mode character object. Please specify here whether you want to
#'   search for colocalizations between miRNAs & miRNAs (="MM") or between
#'   miRNAs & RBPs (="MR") or between miRNAs & a custom position objects
#'   (="MC"). See in cust.input which information this object has to contain.
#' @param cust.input A dataframe containing positional information of any given
#'   object. The dataframe has to have at least 4 columns with the following
#'   information and structure: Ensembl Transcript IDs in the first column (eg.
#'   ENSMUST00000046533), the Motif or Seed specifier of the custom object in
#'   the second column (eg. "Syncrip") and the start and stop position relative
#'   to the 3' UTR start in the third and fourth column respecitively (eg. col3
#'   = 20 & col4 = 27, meaning position 1020 to 1027 if the first NT of the
#'   3'UTR is at position 1001)
#' @param species character object. Can be "human", "mouse" or "rat" at the
#'   moment
#' @param expressed.miRNAs A vector of expressed miRNAs (miRBase IDs), or
#'   specific miRNAs to search for
#' @param co_dist_min Minimum number of nucleotides to be located between two
#'   miRNA-seeds / KD-12mers and taken into consideration for the search. By
#'   specifying co_dist_min and co_dist_max, the user can search for certain
#'   types of functional miRNA colocalization. Mutual miRNA blocking is
#'   suggested to happen within a distance of 0-7Nt. Cooperative repression of
#'   two miRNAs is reported to happen within a distance of 8-40Nt, and at bigger
#'   distance, it is suggested to miRNAs act independently an given targets. If
#'   both, co_dist_min and co_dist_max are left blank, enrichment of miRNA-pairs
#'   on the whole transcript will be reported.
#' @param co_dist_max Maximum number of nucleotides to be located between two
#'   miRNA-seeds / KD-12mers and taken into consideration for the search. For
#'   further info, see co_dist_min.
#' @param BP Pass `BiocParallel::MulticoreParam(ncores, progressbar=TRUE)` to
#'   enable multithreading.
#' @param KD_sites Possibility to use prediced KS-sites for the search. Not yet.
#' @param min_pairs Minimum number of required miRNA-pairs to be displayed in
#'   the result table. Default = 5
#'
#' @return a data.table of all colocalizations including statistics
#'
#' @import GenomicRanges IRanges data.table
#' @export
#'
#' 
#'  
#'    
findCoLoc <- function(expressed.genes=NULL, expressed.miRNAs=NULL, search.mode = c("MM","MR","MC"), cust.input = NULL, species = c("human","mouse","rat"), co_dist_min = NULL, co_dist_max = NULL, 
                      BP= NULL, KD_Sites = FALSE, min_pairs = 5){
  library(GenomicRanges)
  library(IRanges)
  library(data.table)


  if(any(!(search.mode %in% c("MM","MR","MC"))))
    stop("search.mode has to be specified with one of the given input options")
     
  if(!(is.null(cust.input)) && search.mode != "MC")   
     message("The custom input object won't be used in a miRNA-miRNA or miRNA-RBP colocalization search")
  
  if(search.mode == "MC" && is.null(cust.input))
    stop("A custom input object has to be specified if a miRNA-Custom.Object colocalization search is specified in the search.mode")
    
  if(is.null(co_dist_min) & is.null(co_dist_max))
    stop("Either a minium or maximum colocalization distance has to be chosen. If you want to get all pairs, set 'co_dist_min = 0'")
  
  if(co_dist_min >= co_dist_max)
    stop("The specified maximal distance between two motifs has to be bigger than the minimal distance")
  
  if (KD_Sites) 
    stop("not yet")
  else{
    
    #Download Targetscan Positions
    m <- .loadTargetscanPos(species)
    
    #filter for expressed genes
    if(!is.null(expressed.genes)){
      m <- m[m$`Gene ID` %in% expressed.genes,]
    }
    
    #Download miRNA families
    miRFam <- .getMirfamilies(species)
    mirvec <- CharacterList(lapply(split(miRFam$`MiRBase ID`, miRFam$`Seed+m8`),unique))
    
    
    #filter for expressed miRNAs
    if(!is.null(expressed.miRNAs)){
      mirvec <- mirvec[sapply(mirvec, function(y) any(expressed.miRNAs %in% y))]
    }
      
    # include miRnames
    # sollte auch so gehen vapply(mirvec,FUN.VALUE = character(length = 1),FUN=function(x){ paste(x, collapse = ", ") })
    ff <- sapply(mirvec,FUN=function(x){ paste(x, collapse = ", ") })
    ff <- as.data.frame(ff,row.names=names(ff))
    colnames(ff) <- "miRNAs"
    ff <- as.data.table(ff,keep.rownames = TRUE)
    ff <- merge(ff,miRFam[,c("miR family","Seed+m8")], by.x = "rn",by.y ="Seed+m8", all = FALSE)
    colnames(ff)[colnames(ff)=="rn"] <- "Seed+m8"
    ff <- unique(ff,by = "Seed+m8")
      
    #miRNA Position DataFrame preparation with optionally filtered miRs
    m <- merge(m, ff, by.x = "miR Family", by.y = "miR family", all = FALSE)
      

    if(length(m$`Gene ID`) == 0)
      stop("No miRNA pairs found in expressed genes")
  
      
    #Prepare the GRanges Object
    mirs <- GRanges(seqnames=m$`Transcript ID`, IRanges(m$`UTR start`, m$`UTR end`), miRNA=m$`miR Family`, seed=m$`Seed+m8`)
    
    dt <- switch( search.mode,
                      MM = {
                        dt <- .findcoloc1mir(co_dist_min,co_dist_max,mirs)
                        dt},
                      MR = {
                        RBPPos <- .loadRBPPos(species)
                        r <- merge(RBPPos,m[,c("Transcript ID","Gene ID")], by = "Transcript ID")
                        rbps <- GRanges(seqnames=r$`Transcript ID`, IRanges(r$start, r$stop), motif=r$Motif_ID)
                        dt <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs,rbps)                
                        dt},
                      MC = {
                        r <- cust.input
                        colnames(r)[1:4] <- c("Transcript ID","Motif_ID","start","stop")
                        r <- merge(r,m[,c("Transcript ID","Gene ID")], by = "Transcript ID")
                        cust <- GRanges(seqnames=r$`Transcript ID`, IRanges(r$start, r$stop), motif=r$Motif_ID)
                        dt <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs,cust)                
                        dt
                      },
                      stop("This stop actually shouldn't happen")
    )
    
    
    # Applying Statistics
    dt$p.value.enrich <- apply(dt,1,function(dt){
      m <- matrix(as.numeric(dt[c(3, 4, 5, 6)]), ncol = 2)
      f <- fisher.test(as.table(m), alt="greater")
      return(f$p.value)
    })
    dt$FDR <- p.adjust(dt$p.value.enrich, method = "BH")
    
    # Get the names
    m <- m[,c("Seed+m8","miRNAs")]
    m <- unique(m,by = "Seed+m8")
    dt <- merge(dt,m, by.x = "miRNA", by.y = "Seed+m8",all.x = TRUE)
    
    dt <- switch( search.mode,
                  MM = {
                    dt <- merge(dt,m, by.x = "Partner", by.y = "Seed+m8",all.x = TRUE)
                    colnames(dt)[colnames(dt)=="miRNAs.x"] <- "miRNA_mir.names"
                    colnames(dt)[colnames(dt)=="miRNAs.y"] <- "Partner_mir.names"
                    dt},
                  MR = {
                    colnames(dt)[colnames(dt)=="miRNAs"] <- "miRNA_names"
                    RBPPos <- RBPPos[,c("Motif_ID","RBP_mus_all","RBP_mus_direct")]
                    RBPPos <- unique(RBPPos,by = "Motif_ID")
                    dt <- merge(dt,RBPPos,by.x = "Partner",by.y = "Motif_ID", all.x = TRUE)
                    colnames(dt)[colnames(dt)=="RBP_mus_direct"] <- "RBP_names.direct"
                    colnames(dt)[colnames(dt)=="RBP_mus_all"] <- "RBP_names.all"
                    setcolorder(dt,c(1:7,9,8))
                    dt},
                  MC = {
                    colnames(dt)[colnames(dt)=="miRNAs"] <- "miRNA_names"
                    colnames(dt)[colnames(dt)=="Motif_ID"] <- "Custom.Motif"
                    dt
                  },
                  stop("This stop actually shouldn't happen")
    )
    
    dt <- dt[dt$overlap >= min_pairs,]
    dt <- setcolorder(dt,"miRNA")
    dt[order(p.value.enrich),]  
    
  }
}    



#' Find miRNA colocalizations in a subset of genes
#'
#' Finds pairs of miRNAs that are located on the same gene or transcript within
#' a distance that the user can specify. Fisher's enrichment test is applied in
#' order to determine which miRNA-pairs are significantly enriched in the
#' particular universe of expressed genes.
#'
#' @param genes_in_set A vector of genes, that are a subset of the expressed
#'   genes. If specified, the colocalization function will test for an
#'   enrichment of colocalizations in these "genes_in_set" with regards to the
#'   gene universe ("expressed.genes").
#' @param expressed.genes A vector of expressed genes (at the moment EnsemblID,
#'   maybe give option?) as a universe to search in (including the genes_in_set)
#' @param search.mode character object. Please specify here whether you want to
#'   search for colocalizations between miRNAs & miRNAs (="MM") or between
#'   miRNAs & RBPs (="MR") or between miRNAs & a custom position objects
#'   (="MC"). See in cust.input which information this object has to contain.
#' @param cust.input A dataframe containing positional information of any given
#'   object. The dataframe has to have at least 4 columns with the following
#'   information and structure: Ensembl Transcript IDs in the first column (eg.
#'   ENSMUST00000046533), the Motif or Seed specifier of the custom object in
#'   the second column (eg. "Syncrip") and the start and stop position relative
#'   to the 3' UTR start in the third and fourth column respecitively (eg. col3
#'   = 20 & col4 = 27, meaning position 1020 to 1027 if the first NT of the
#'   3'UTR is at position 1001)
#' @param species character object. Can be "human", "mouse" or "rat" at the
#'   moment
#' @param expressed.miRNAs A vector of expressed miRNAs (miRBase IDs)
#' @param co_dist_min Minimum number of nucleotides to be located between two
#'   miRNA-seeds / KD-12mers and taken into consideration for the search. By
#'   specifying co_dist_min and co_dist_max, the user can search for certain
#'   types of functional miRNA colocalization. Mutual miRNA blocking is
#'   suggested to happen within a distance of 0-7Nt. Cooperative repression of
#'   two miRNAs is reported to happen within a distance of 8-40Nt, and at bigger
#'   distance, it is suggested to miRNAs act independently an given targets. If
#'   both, co_dist_min and co_dist_max are left blank, enrichment of miRNA-pairs
#'   on the whole transcript will be reported.
#' @param co_dist_max Maximum number of nucleotides to be located between two
#'   miRNA-seeds / KD-12mers and taken into consideration for the search. For
#'   further info, see co_dist_min.
#' @param Possibility to use prediced KS-sites for the search. Not yet.
#' @param min_pairs_sub Minimum number of required miRNA-Partner-pairs in the
#'   subset of genes to be displayed in the result table. Default = 5
#'
#' @return a data.table of all colocalizations including statistics
#'
#' @import GenomicRanges IRanges data.table
#' @export
#'
#'
#'
#' 
#' 
#' 
#' 
findCoLocsubset <- function(genes_in_set=NULL, expressed.genes=NULL, search.mode = c("MM","MR","MC"), cust.input = NULL, species = c("human","mouse","rat"), expressed.miRNAs=NULL, co_dist_min = NULL, co_dist_max = NULL, KD_Sites = FALSE, min_pairs_sub = 5){
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  
  if(any(!(search.mode %in% c("MM","MR","MC"))))
    stop("search.mode has to be specified with one of the given input parameters")
  
  if(!(is.null(cust.input)) && search.mode != "MC")   
    message("The custom input object won't be used in a miRNA-miRNA or miRNA-RBP colocalization search")
  
  if(search.mode == "MC" && is.null(cust.input))
    stop("A custom input object has to be specified if a miRNA-Custom.Object colocalization search is specified in the search.mode")
  
  if(is.null(co_dist_min) & is.null(co_dist_max))
    stop("Either a minium or maximum colocalization distance has to be chosen. If you want to get all pairs, set 'co_dist_min = 0'")
  
  if(any(!(genes_in_set %in% expressed.genes)))
    stop("All specified subset genes have to be present in the gene universe")
  
  
  
  if (KD_Sites) message("not yet")
  else{
    
    #Download Targetscan Positions
    m <- .loadTargetscanPos(species)
    
    #filter for expressed genes
    if(!is.null(expressed.genes)){
      m <- m[m$`Gene ID` %in% expressed.genes,]
    }
    
    #Download miRNA families
    miRFam <- .getMirfamilies(species)
    mirvec <- CharacterList(lapply(split(miRFam$`MiRBase ID`, miRFam$`Seed+m8`),unique))
    
    
    #filter for expressed miRNAs
    ## is that correct?
    if(!is.null(expressed.miRNAs)){
      mirvec <- mirvec[sapply(mirvec, function(y) any(expressed.miRNAs %in% y))]
    }
    
    # include miRnames
    ff <- sapply(mirvec,FUN=function(x){ paste(x, collapse = ", ") })
    ff <- as.data.frame(ff,row.names=names(ff))
    colnames(ff) <- "miRNAs"
    ff <- as.data.table(ff,keep.rownames = TRUE)
    ff <- merge(ff,miRFam[,c("miR family","Seed+m8")], by.x = "rn",by.y ="Seed+m8", all.x = TRUE)
    colnames(ff)[colnames(ff)=="rn"] <- "Seed+m8"
    ff <- unique(ff,by = "Seed+m8")
    
    #miRNA Position DataFrame preparation with optionally filtered miRs
    m <- merge(m, ff, by.x = "miR Family", by.y = "miR family", all = FALSE)
    
    
    if(length(m$`Gene ID`) == 0)
      stop("No miRNA pairs found in expressed genes")
    
    
    ## Possibiliy would be the following:
    m <- as.data.table(m)
    m <- m[m$`Gene ID` %in% genes_in_set,"gene_info":= "in_set"]
    m <- m[!(m$`Gene ID` %in% genes_in_set),"gene_info":= "back"]
    
    #Prepare the GRanges Object
    mirs <- GRanges(seqnames=m$`Transcript ID`, IRanges(m$`UTR start`, m$`UTR end`), miRNA=m$`miR Family`, seed=m$`Seed+m8`, info = m$gene_info)
    mirs1 <- mirs[mirs$info == "in_set",]
    mirs2 <- mirs[mirs$info == "back",]
    
    ###### 
    
    dt <- switch( search.mode,
                  MM = {
                    dt1 <- .findcoloc1mir(co_dist_min,co_dist_max,mirs1)
                    colnames(dt1)[-c(1,2)] <- paste(colnames(dt1)[-c(1,2)], "in.subset", sep = "_")
                    
                    dt2 <- .findcoloc1mir(co_dist_min,co_dist_max,mirs2)
                    dt <- merge(dt1,dt2,by = c("miRNA","Partner"), all = TRUE)
                    dt[is.na(dt)] <- 0
                    dt},
                  MR = {
                    RBPPos <- .loadRBPPos(species)
                    r <- merge(RBPPos,m[,c("Transcript ID","Gene ID", "gene_info")], by = "Transcript ID")
                    rbps <- GRanges(seqnames=r$`Transcript ID`, IRanges(r$start, r$stop), motif=r$Motif_ID, info = r$gene_info)
                    
                    rbps1 <- rbps[rbps$info == "in_set",]
                    rbps2 <- rbps[rbps$info == "back" ,]
                    
                    dt1 <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs1,rbps1)     
                    colnames(dt1)[-c(1,2)] <- paste(colnames(dt1)[-c(1,2)], "in.subset", sep = "_")
                    
                    dt2 <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs2, rbps2)
                    
                    dt <- merge(dt1,dt2,by = c("miRNA","Partner"), all = TRUE)
                    dt[is.na(dt)] <- 0
                    dt},
                  MC = {
                    r <- cust.input
                    colnames(r)[1:4] <- c("Transcript ID","Motif_ID","start","stop")
                    r <- merge(r,m[,c("Transcript ID","Gene ID", "gene_info")], by = "Transcript ID")
                    cust <- GRanges(seqnames=r$`Transcript ID`, IRanges(r$start, r$stop), motif=r$Motif_ID, info = r$gene_info )
                    
                    cust1 <- cust[cust$info == "in_set",]
                    cust2 <- cust[cust$info == "back" ,]
                    
                    dt1 <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs1,cust1)     
                    colnames(dt1)[-c(1,2)] <- paste(colnames(dt1)[-c(1,2)], "in.subset", sep = "_")
                    
                    dt2 <- .findcoloc1rbp(co_dist_min,co_dist_max,mirs2, cust2)
                    
                    dt <- merge(dt1,dt2,by = c("miRNA","Partner"), all = TRUE)
                    dt[is.na(dt)] <- 0
                    dt},
                  stop("This stop actually shouldn't happen")
    )
    
    ####
    # We build three different statistics:
    # 1) Check for the enrichment of a specific pair in the subset vs the background.
    # 2) Check whether the specified microRNA has more colocolizations with this partner than
    #    with any other.
    # 3) Check whether the specified Partner has more colocolizations with this microRNA than
    #    with any other.
    
    
    # 1)  
    dt$overlap.enrich <- apply(dt,1,function(dt){
      m <- matrix(as.numeric(dt[c(3, 6, 7, 10)]), ncol = 2)
      f <- fisher.test(as.table(m), alt="greater")
      return(f$p.value)
    })
    dt$FDR.overlap.enrich <- p.adjust(dt$overlap.enrich, method = "BH")
    
    
    # 2)   
    dt$miRNA_pair.enrich <- apply(dt,1,function(dt){
      m <- matrix(as.numeric(dt[c(3, 4, 7, 8)]), ncol = 2)
      f <- fisher.test(as.table(m), alt="greater")
      return(f$p.value)
    })
    dt$FDR.miRNA_pair.enrichr <- p.adjust(dt$miRNA_pair.enrich, method = "BH")
    
    
    # 3)
    dt$Partner_pair.enrich <- apply(dt,1,function(dt){
      m <- matrix(as.numeric(dt[c(3, 5, 7, 9)]), ncol = 2)
      f <- fisher.test(as.table(m), alt="greater")
      return(f$p.value)
    })
    dt$FDR.Partner_pair.enrichr <- p.adjust(dt$Partner_pair.enrich, method = "BH")
    
    
    # Get the names
    m <- m[,c("Seed+m8","miRNAs")]
    m <- unique(m,by = "Seed+m8")
    dt <- merge(dt,m, by.x = "miRNA", by.y = "Seed+m8",all.x = TRUE)
    
    dt <- switch( search.mode,
                  MM = {
                    dt <- merge(dt,m, by.x = "Partner", by.y = "Seed+m8",all.x = TRUE)
                    colnames(dt)[colnames(dt)=="miRNAs.x"] <- "miRNA_mir.names"
                    colnames(dt)[colnames(dt)=="miRNAs.y"] <- "Partner_mir.names"
                    dt},
                  MR = {
                    colnames(dt)[colnames(dt)=="miRNAs"] <- "miRNA_names"
                    RBPPos <- RBPPos[,c("Motif_ID","RBP_mus_all","RBP_mus_direct")]
                    RBPPos <- unique(RBPPos,by = "Motif_ID")
                    dt <- merge(dt,RBPPos,by.x = "Partner",by.y = "Motif_ID", all.x = TRUE)
                    colnames(dt)[colnames(dt)=="RBP_mus_direct"] <- "RBP_names.direct"
                    colnames(dt)[colnames(dt)=="RBP_mus_all"] <- "RBP_names.all"
                    setcolorder(dt,c(1:17,19,18))
                    dt},
                  MC = {
                    colnames(dt)[colnames(dt)=="miRNAs"] <- "miRNA_names"
                    colnames(dt)[colnames(dt)=="Motif_ID"] <- "Custom.Motif"
                    dt
                  },
                  stop("This stop actually shouldn't happen")
    )
    
    dt <- dt[dt$overlap_in.subset >= min_pairs_sub,]
    dt <- setcolorder(dt,"miRNA")
    dt[order(overlap.enrich),]  
    
  }
}    







.findcoloc1mir <- function(min.dist,max.dist,Pos_object) {
  df <- sapply( split(Pos_object, Pos_object$seed), FUN=function(mirs_all){
    # each mirs contains all the binding sites of a single miRNA
    # this means basically that a second argument (the mirs_all) is passed to the lapply function
    ll <- lapply(split(Pos_object, Pos_object$seed), mi=mirs_all, FUN=function(x,mi){
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
  dt
}









.findcoloc1rbp<- function(min.dist,max.dist,mirs,Pos_Object2) {
 
  # This should be miRNA centered, correct? Would that matter?
  
  df <- sapply( split(mirs, mirs$seed), FUN=function(mirs_all){
    # mirs contains all the binding sites of a single miRNA
    ll <- lapply(split(rbps, rbps$motif), mi=mirs_all, FUN=function(x,mi){
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
  dt

}












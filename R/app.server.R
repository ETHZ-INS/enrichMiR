#' enrichMiR.server
#'
#' @param modlists A named list of `CompressedKdModelList`
#' @param targetlists An optional list of aggregated targets, named as `modlists`
#' @param ensdbs A named list of `ensembldb` objects (with names corresponding to those 
#' of `modlists`)
#' @param genomes A named list of `BSgenome` objects (with names corresponding to those 
#' of `modlists`)
#' @param maxCacheSize Maximum cache size in bytes.
#'
#' @return A shiny server function
#' @importFrom digest digest
#' @export
enrichMiR.server <- function(modlists, targetlists=list(), ensdbs=list(), genomes=list(),
                             maxCacheSize=100*10^6 ){
  library(DT)
  library(stringr)
  library(Biostrings)
  library(ggplot2)
  
  dtwrapper <- function(d, pageLength=25){
    datatable( d, filter="top", class="compact", extensions=c("Buttons","ColReorder"),
               options=list(pageLength=pageLength, dom = "fltBip", rownames=FALSE,
                            colReorder=TRUE, 
                            buttons=c('copy', 'csv', 'excel', 'csvHtml5') ) )
  }
  
  
  function(input, output, session){
    
    ##############################
    ## initialize inputs
    
    updateSelectizeInput(session, "mirlist", choices=names(modlists))
    updateSelectizeInput(session, "annotation", choices=names(ensdbs))
    
    observe({
      # when the choice of collection changes, update the annotation to
      # use the same genome
      if(!is.null(input$mirlist) && input$mirlist!="")
        updateSelectizeInput(session, "annotation", selected=input$mirlist)
    })
    
    
    ##############################
    ## select collection
    
    allmods <- reactive({ # all models from collection
      if(is.null(input$mirlist)) return(NULL)
      modlists[[input$mirlist]]
    })
    
    # prints a summary of the model collection
    output$collection_summary <- renderPrint({
      if(is.null(allmods())) return(NULL)
      summary(allmods())
    })
    
    observe({ # when the selected collection changes, update the miRNA selection inputs
      updateSelectizeInput(session, "mirnas", choices=names(allmods()), server=TRUE)
      updateSelectizeInput(session, "mirna", choices=names(allmods()), server=TRUE)
    })
    
    ##############################    
    ## scan specific sequence
    
    ## transcript selection
    
    sel_ensdb <- reactive({ # the ensembldb for the selected genome
      if(is.null(input$annotation) || input$annotation=="" || 
         !(input$annotation %in% names(ensdbs))) return(NULL)
      ensdbs[[input$annotation]]
    })
    
    allgenes <- reactive({ # all genes in the selected genome
      if(is.null(sel_ensdb())) return(NULL)
      g <- genes( sel_ensdb(), columns="gene_name",
                  return.type="data.frame")
      paste(g[,1], g[,2])
    })
    
    selgene <- reactive({ # selected gene id
      if(is.null(input$gene) || input$gene=="") return(NULL)
      strsplit(input$gene," ")[[1]][2]
    })
    
    alltxs <- reactive({ # all tx from selected gene
      if(is.null(selgene())) return(NULL)
      tx <- transcripts(sel_ensdb(), columns=c("tx_id","tx_biotype"),
                        filter=~gene_id==selgene(), return.type="data.frame")
      paste0(tx$tx_id, " (", tx$tx_biotype,")")
    })
    
    seltx <- reactive({ # the selected transcript
      if(is.null(selgene()) || is.null(input$transcript) || 
         input$transcript=="") return(NULL)
      tx <- strsplit(input$transcript, " ")[[1]][[1]]
      if(tx=="" || is.na(tx)) return(NULL)
      tx
    })
    
    # when the ensembldb is updated, update the gene input
    observe(updateSelectizeInput(session, "gene", choices=allgenes(), server=TRUE))
    # when the gene selection is updated, update the transcript input
    observe(updateSelectizeInput(session, "transcript", choices=alltxs()))
    
    # takes a genome package name as input, and returns the genome
    getGenome <- function(x){
      if(is.character(x)){
        library(x, character.only=TRUE)
        x <- get(x)
      }
      seqlevels(x) <- gsub("^chr","",seqlevels(x))
      x
    }
    
    seqs <- reactive({ # returns the selected sequence(s)
      if(is.null(selgene())) return(NULL)
      if(is.null(seltx())){
        gid <- selgene()
        filt <- ~gene_id==gid
      }else{
        txid <- seltx()
        filt <- ~tx_id==txid
      }
      db <- sel_ensdb()
      if(input$utr_only){
        gr <- suppressWarnings(threeUTRsByTranscript(db, filter=filt))
      }else{
        gr <- exonsBy(db, by="tx", filter=filt)
      }
      if(length(gr)==0) return(NULL)
      if(is.null(genomes[[input$annotation]])) return(NULL)
      seqs <- extractTranscriptSeqs(getGenome(genomes[[input$annotation]]), gr)
      seqs <- seqs[lengths(seqs)>6]
      if(length(seqs)==0) return(NULL)
      seqs
    })
    
    output$tx_overview <- renderPrint({ # overview of the selected transcript
      if(is.null(seqs())) return(NULL)
      if(length(seqs())==0 || length(seqs())>1) return(seqs())
      return(seqs()[[1]])
    })
    
    ## end transcript selection

    ## custom sequence
    
    customTarget <- reactive({
      if(is.null(input$customseq) || input$customseq=="") return(NULL)
      seqtype <- suppressWarnings(enrichMiR:::.guessSeqType(input$customseq))
      seq <- input$customseq
      if(input$circular) seq <- paste0(seq,substr(seq,1,min(nchar(seq),11)))
      if(seqtype=="DNA") return(DNAString(seq))
      return(RNAString(seq))
    })
    
    output$custom_info <- renderPrint({ # overview of the custom sequence
      if(is.null(input$customseq)) return("")
      out <- capture.output(customTarget())
      if(input$circular) out <- c("Circularized sequence: the first 11nt are pasted to ",
                                  "the end of the sequence.\n", out)
      cat(out)
    })

    target <- reactive({ # target subject sequence
      if(input$subjet_type=="custom"){
        return(as.character(DNAStringSet(customTarget())))
      }   
      if(is.null(seqs()) || length(seqs())>1) return(NULL)
      as.character(seqs())
    })
    
    observeEvent(input$rndseq, { # generate random sequence
      updateTextAreaInput(session, "customseq",
                          value=paste(sample(c("A","C","G","T"), size = 3000, replace=TRUE), collapse=""))
    })
    
    ## Select miRNAs for scanning
    
    observeEvent(input$mirnas_confident, {
      if(is.null(allmods())) return(NULL)
      cons <- conservation(allmods())
      if(all(is.na(cons))) return(NULL)
      updateSelectizeInput(session, "mirnas", selected=names(cons)[as.numeric(cons)>1])
    })
    observeEvent(input$mirnas_mammals, {
      if(is.null(allmods())) return(NULL)
      cons <- conservation(allmods())
      if(all(is.na(cons))) return(NULL)
      updateSelectizeInput(session, "mirnas", selected=names(cons)[as.numeric(cons)>2])
    })
    observeEvent(input$mirnas_vert, {
      if(is.null(allmods())) return(NULL)
      cons <- conservation(allmods())
      if(all(is.na(cons))) return(NULL)
      updateSelectizeInput(session, "mirnas", selected=names(cons)[as.numeric(cons)>3])
    })
    observeEvent(input$mirnas_clear, {
      updateSelectizeInput(session, "mirnas", selected="")
    })
    
    selmods <- reactive({ # models selected for scanning
      if(is.null(allmods())) return(NULL)
      if(input$mirnas_all) return(allmods())
      if(is.null(input$mirnas)) return(NULL)
      allmods()[input$mirnas]
    })
    
    ## Begin scan and results caching
    
    cached.hits <- reactiveValues() # actual and past scanning results are stored in this object
    
    cached.checksums <- reactive({
      ch <- reactiveValuesToList(cached.hits)
      ch <- ch[!sapply(ch, is.null)]
      ch <- lapply(ch, FUN=function(x){
        x[c("target","size","time","last","nsel","sel")]
      })
    })
    current.cs <- reactiveVal()
    
    hits <- reactive({ # the results currently loaded are stored in this object
      if(is.null(current.cs())) return(NULL)
      if(current.cs() %in% names(cached.checksums())) return(cached.hits[[current.cs()]])
      NULL
    })
    
    checksum <- reactive({ # generate a unique hash for the given input
      paste( digest::digest(selmods()),
             digest::digest(list(target=target(), shadow=input$shadow,
                                 keepMatchSeq=input$keepMatchSeq, 
                                 minDist=input$minDist,
                                 scanNonCanonical=input$scanNonCanonical))
      )
    })
    
    cache.size <- reactive({
      ch <- cached.checksums()
      if(is.null(ch) || length(ch)==0) return(0)
      sum(sapply(ch, FUN=function(x) as.numeric(x$size)))
    })
    
    cleanCache <- function(){ # remove last-used results when over the cache size limit
      cs <- isolate(cached.checksums())
      if(length(cs)<3 || as.numeric(cache.size())<maxCacheSize) return(NULL)
      cs <- cs[order(sapply(cs, FUN=function(x) x$last), decreasing=TRUE)]
      sizes <- sapply(cs, FUN=function(x) as.numeric(x$size))
      while(length(cs)>2 & sum(sizes)>maxCacheSize){
        cached.hits[[rev(names(cs))[1]]] <- NULL
        cs <- cs[-length(cs)]
        sizes <- sizes[-length(sizes)]
      }
    }
    
    observeEvent(input$scan, { # actual scanning
      if(is.null(selmods()) || is.null(target()) || nchar(target())==0) 
        return(NULL)
      cs <- checksum()
      if(!(cs %in% names(cached.checksums()))){
        msg <- paste0("Scanning sequence for ",length(selmods())," miRNAS")
        detail <- NULL
        if(length(selmods())>4) detail <- "This might take a while..."
        if(input$circular) detail <- "'Ribosomal Shadow' is ignored when scanning circRNAs"
        withProgress(message=msg, detail=detail, value=1, max=3, {
          cached.hits[[cs]]$hits <- findSeedMatches( target(), selmods(),
                                                     keepMatchSeq=input$keepmatchseq,
                                                     minDist=input$minDist,
                                                     shadow=ifelse(input$circular,0,input$shadow),
                                                     max.noncanonical.motifs=ifelse(input$scanNonCanonical,Inf,0),
                                                     BP=MulticoreParam(2, progressbar=TRUE) )
          if(length(cached.hits[[cs]])>0){
            cached.hits[[cs]]$hits$miRNA <- cached.hits[[cs]]$hits$seed
            cached.hits[[cs]]$hits$seed <- NULL
          }
        })
        cached.hits[[cs]]$cs <- cs
        cached.hits[[cs]]$last <- cached.hits[[cs]]$time <- Sys.time()
        cached.hits[[cs]]$size <- object.size(cached.hits[[cs]]$hits)
        cached.hits[[cs]]$nsel <- nm <- length(selmods())
        cached.hits[[cs]]$sel <- ifelse(nm>1,paste(nm,"models"),input$mirnas)
        cached.hits[[cs]]$target_length <- nchar(target())
        if(input$subjet_type=="custom"){
          cached.hits[[cs]]$target <- "custom sequence"
        }else{
          cached.hits[[cs]]$target <- paste(input$gene, "-", seltx(),
                                            ifelse(input$utr_only, "(UTR)",""))
        }
        cleanCache()
      }
      current.cs(cs)
    })
    
    output$scan_target <- renderText({
      if(is.null(current.cs()) || is.null(cached.hits[[current.cs()]])) return(NULL)
      paste("Scan results in: ", cached.hits[[current.cs()]]$target)
    })
    
    output$cache.info <- renderText({
      if(cache.size()==0) return("Cache empty.")
      paste0(length(cached.checksums()), " results cached (",
             round(cache.size()/1024^2,3)," Mb)")
    })
    
    output$cached.results <- renderUI({
      ch <- cached.checksums()
      ch2 <- names(ch)
      names(ch2) <- sapply(ch, FUN=function(x){
        paste0(x$time, ": ", x$sel, " on ", x$target, " (", format(x$size,units="Kb"),")")
      })
      radioButtons("selected.cache", "Cached results", choices=ch2)
    })
    
    observeEvent(input$loadCache, {
      if(is.null(input$selected.cache)) return(NULL)
      current.cs(input$selected.cache)
    })
    
    observeEvent(input$deleteCache, {
      if(is.null(input$selected.cache)) return(NULL)
      cached.hits[[input$selected.cache]] <- NULL
      if(current.cs()==input$selected.cache) current.cs(NULL)
    })
    
    output$hits_table <- renderDT({ # prints the current hits
      if(is.null(hits()$hits)) return(NULL)
      h <- as.data.frame(hits()$hits)
      h <- h[,setdiff(colnames(h), c("seqnames","width","strand") )]
      if(hits()$nsel == 1) h$miRNA <- NULL
      dtwrapper(h)
    })
    
    ## end scan hits and cache 
    
    output$manhattan <- renderPlotly({
      if(is.null(hits()$hits)) return(NULL)
      h <- as.data.frame(sort(hits()$hits))
      tt <- sort(table(h$miRNA), decreasing=TRUE)
      mirs <- names(tt)[seq_len(min(input$manhattan_n, length(tt)))]
      h <- h[h$miRNA %in% mirs,,drop=FALSE]
      if(!is.null(input$manhattan_ordinal) && input$manhattan_ordinal){
        h$position <- seq_len(nrow(h))
        xlab <- "Position (ordinal)"
        xlim <- c(0,nchar(hits()$target_length))
      }else{
        h$position <- round(rowMeans(h[,2:3]))
        xlab <- "Position (nt) in sequence"
        xlim <- NULL
      }
      if("sequence" %in% colnames(h)){
        p <- ggplot(h, aes(position, -log_kd, colour=miRNA, seq=sequence, kmer_type=type))
      }else{
        p <- ggplot(h, aes(position, -log_kd, colour=miRNA, kmer_type=type))
      }
      p <- p + geom_hline(yintercept=1.5, linetype="dashed", color = "red", size=1) + 
        geom_point(size=2) + xlab(xlab) + expand_limits(x=xlim, y=0)
      ggplotly(p)
    })
    
    ##############################
    ## miRNA-centric tab
    
    mod <- reactive({ # the currently-selected KdModel
      if(is.null(allmods()) || is.null(input$mirna)) return(NULL)
      allmods()[[input$mirna]]
    })
    
    output$modconservation <- renderText({
      if(is.null(mod())) return(NULL)
      as.character(conservation(mod()))
    })
    
    output$mirbase_link <- renderUI({
      tags$a(href=paste0("http://www.mirbase.org/textsearch.shtml?q=", input$mirna), 
             "miRBase", target="_blank")
    })
    
    output$modplot <- renderPlot({ # affinity plot
      if(is.null(mod())) return(NULL)
      plotKdModel(mod())
    }, height=reactive(input$modplot_height))
    
    txs <- reactive({ # the tx to gene symbol table for the current annotation
      if(is.null(input$mirlist) || input$mirlist=="" || 
         !(input$mirlist %in% names(ensdbs))) return(NULL)
      db <- ensdbs[[input$mirlist]]
      if(is.null(db)) return(NULL)
      tx <- mcols(transcripts(db, c("tx_id","gene_id","tx_biotype")))
      tx <- merge(tx,mcols(genes(db, c("gene_id","symbol"))), by="gene_id")
      as.data.frame(tx[,c("symbol","tx_id","tx_biotype")])
    })
    
    output$mirna_targets <- renderDT({
      tl <- paste0(input$mirlist, ifelse(input$targetlist_utronly, ".utrs", ".full"))
      if(is.null(targetlists[[tl]]) || is.null(mod())) return(NULL)
      d <- targetlists[[tl]][[input$mirna]]
      d$log_kd <- abs(d$log_kd/100)
      d$log_kd.canonical <- abs(d$log_kd.canonical/100)
      if(!is.null(txs())){
        d <- merge(txs(), d, by.x="tx_id", by.y="transcript")
        if(input$targetlist_gene)
          d <- aggregate(d[,grep("mer|log_kd",colnames(d))], d[,c("symbol","seed")], na.rm=TRUE, FUN=max)
      }
      d$log_kd <- -d$log_kd
      d$log_kd.canonical <- -d$log_kd.canonical
      d <- d[which(d$log_kd <= -1.5),]
      colnames(d) <- gsub("^n\\.","",colnames(d))
      dtwrapper(d)
    })
    
  }
}

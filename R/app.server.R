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
#' @export
enrichMiR.server <- function(modlists, targetlists=list(), ensdbs=list(), genomes=list(),
                             maxCacheSize=100*10^6){
  
  function(input, output, session){
    
    ##############################
    ## initialize inputs
    
    updateSelectizeInput(session, "mirlist", choices=names(modlists))
    updateSelectizeInput(session, "annotation", choices=names(ensdbs))
    
    
    ##############################
    ## select collection
    
    allmods <- reactive({ # all models from collection
      if(is.null(input$mirlist)) return(NULL)
      modlists[[input$mirlist]]
    })
    
    output$collection_summary <- renderPrint({
      if(is.null(allmods())) return(NULL)
      summary(allmods())
    })
    
    observe({
      updateSelectizeInput(session, "mirnas", choices=names(allmods()), server=TRUE)
      updateSelectizeInput(session, "mirna", choices=names(allmods()), server=TRUE)
    })
        
    ##############################    
    ## scan specific sequence
    
    ## transcript selection
    
    ensdb <- reactive({ # active ensemblDB
      if(is.null(input$annotation)) return(NULL)
      ensdbs[[input$annotation]]
    })
    
    allgenes <- reactive({
      if(is.null(ensdb())) return(NULL)
      g <- genes( ensdb(), columns="gene_name",
                  return.type="data.frame")
      paste(g[,1], g[,2])
    })
    
    selgene <- reactive({ # selected gene id
      if(is.null(ensdb()) || is.null(input$gene) ||
         input$gene=="") return(NULL)
      strsplit(input$gene," ")[[1]][2]
    })
    
    alltxs <- reactive({ # all tx from selected gene
      if(is.null(selgene())) return(NULL)
      tx <- transcripts(ensdb(),
                        columns=c("tx_id","tx_biotype"),
                        filter=~gene_id==selgene(), return.type="data.frame")
      paste0(tx$tx_id, " (", tx$tx_biotype,")")
    })
    
    seltx <- reactive({
      if(is.null(selgene()) || is.null(input$transcript) || 
         input$transcript=="") return(NULL)
      tx <- strsplit(input$transcript, " ")[[1]][[1]]
      if(tx=="" || is.na(tx)) return(NULL)
      tx
    })
    
    observe(
      updateSelectizeInput(session, "gene", choices=allgenes(), server=TRUE)
    )
    observe(updateSelectizeInput(session, "transcript", choices=alltxs()))
    
    seqs <- reactive({
      if(is.null(selgene())) return(NULL)
      if(is.null(seltx())){
        gid <- selgene()
        filt <- ~gene_id==gid
      }else{
        txid <- seltx()
        filt <- ~tx_id==txid
      }
      if(input$utr_only){
        gr <- suppressWarnings(threeUTRsByTranscript(ensdb(), filter=filt))
      }else{
        gr <- exonsBy(ensdb(), by="tx", filter=filt)
      }
      if(length(gr)==0) return(NULL)
      seqs <- extractTranscriptSeqs(genomes[[input$annotation]], gr)
      seqs <- seqs[lengths(seqs)>6]
      if(length(seqs)==0) return(NULL)
      seqs
    })
    
    output$tx_overview <- renderPrint({
      if(is.null(seqs())) return(NULL)
      if(length(seqs())==0 || length(seqs())>1) return(seqs())
      return(seqs()[[1]])
    })
    
    ## end transcript selection
        
    target <- reactive({ # target subject sequence
      if(input$subjet_type=="custom") return(input$customseq)
      if(is.null(seqs()) || length(seqs())>1) return(NULL)
      as.character(seqs())
    })
    
    observeEvent(input$rndseq, {
      updateTextAreaInput(session, "customseq",
        value=paste(sample(c("A","C","G","T"), size = 3000, replace=TRUE), collapse=""))
    })
    
    selmods <- reactive({ # models selected for scanning
      if(is.null(allmods())) return(NULL)
      if(input$mirnas_all) return(allmods())
      if(is.null(input$mirnas)) return(NULL)
      allmods()[input$mirnas]
    })
    
    ## Begin scan and results caching
    
    cached.hits <- reactiveValues()
    cached.checksums <- reactive({
      ch <- reactiveValuesToList(cached.hits)
      ch <- ch[!sapply(ch, is.null)]
      ch <- lapply(ch, FUN=function(x){
        x[c("target","size","time","last","nsel","sel")]
      })
    })
    current.cs <- reactiveVal()
    
    hits <- reactive({
      if(is.null(current.cs())) return(NULL)
      if(current.cs() %in% names(cached.checksums())) return(cached.hits[[current.cs()]])
      NULL
    })
    
    checksum <- reactive({
      paste( R.cache::getChecksum(selmods()),
             R.cache::getChecksum(list(target=target(), shadow=input$shadow,
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
    
    cleanCache <- function(){
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
        withProgress(message=msg, detail=detail, value=1, max=3, {
          cached.hits[[cs]]$hits <- findSeedMatches( target(), selmods(),
                                        keepMatchSeq=input$keepmatchseq,
                                        minDist=input$minDist,
                                        shadow=input$shadow,
                                        max.noncanonical.motifs=ifelse(input$scanNonCanonical,Inf,0),
                                        BP=SerialParam(progressbar=TRUE) )
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
    
    output$hits_table <- renderDT({
      if(is.null(hits()$hits)) return(NULL)
      h <- as.data.frame(hits()$hits)
      h <- h[,setdiff(colnames(h), c("seqnames","width","strand") )]
      if(hits()$nsel == 1) h$miRNA <- NULL
      datatable( h, options=list( pageLength=25, dom = "fltBip",
                                  buttons=c("csvHtml5","excelHtml5") ),
                 filter="top",
                 extensions=c("Buttons") )
    })
    
    ## end scan hits and cache 
    
    output$manhattan <- renderPlot({
      if(is.null(hits()$hits)) return(NULL)
      h <- as.data.frame(sort(hits()$hits))
      if(input$manhattan_ordinal){
        h$position <- seq_len(nrow(h))
        xlab <- "Position (ordinal)"
      }else{
        h$position <- round(rowMeans(h[,2:3]))
        xlab <- "Position (nt)"
      }
      ggplot(h, aes(position, -log_kd, colour=miRNA)) + geom_point(size=2) + 
        xlab(xlab)
    })
    
    ##############################
    ## miRNA-centric tab

    mod <- reactive({
      if(is.null(allmods()) || is.null(input$mirna)) return(NULL)
      allmods()[[input$mirna]]
    })
    
    output$modplot <- renderPlot({
      if(is.null(mod())) return(NULL)
      plotKdModel(mod())
    }, height=reactive(input$modplot_height))
    
    output$mirna_targets <- renderDT({
      if(is.null(targetlists[[input$mirlist]]) || is.null(mod())) return(NULL)
      d <- targetlists[[input$mirlist]][[input$mirna]]
      d$log_kd <- d$log_kd/100
      d$log_kd.canonical <- d$log_kd.canonical/100
      if(!is.null(ensdbs[[input$mirlist]])){
        tx <- mcols(transcripts(ensdbs[[input$mirlist]], c("tx_id","gene_id","tx_biotype")))
        tx <- merge(tx,mcols(genes(ensdbs[[input$mirlist]], c("gene_id","symbol"))),by="gene_id")
        # tx <- as.data.frame(tx[,c("symbol","gene_id","tx_id","tx_biotype")])
        tx <- as.data.frame(tx[,c("symbol","tx_id")])
        d <- merge(tx, d, by.x="tx_id", by.y="transcript")
      }
      colnames(d) <- gsub("^n\\.","",colnames(d))
      datatable( d, options=list( pageLength=25, dom = "fltBip", rownames=FALSE,
                                  buttons=c("csvHtml5","excelHtml5") ),
                 filter="top", extensions=c("Buttons","ColReorder") )
    })

  }
}
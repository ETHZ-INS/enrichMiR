enrichMiR.server <- function(modlists, targetlists=list(), ensdbs=list(), genomes=list()){
  
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
      if(is.null(allmods()) || is.null(input$mirnas)) return(NULL)
      if(input$mirnas_all) return(allmods())
      allmods()[input$mirnas]
    })
    
    hits <- reactiveValues()
    
    observeEvent(input$scan, { # actual scanning
      if(is.null(selmods()) || is.null(target()) || nchar(target())==0) 
        return(NULL)
      hits$hits <- hits$target <- NULL
      msg <- paste0("Scanning sequence for ",length(selmods())," miRNAS")
      detail <- NULL
      hits$nsel <- length(selmods())
      if(length(selmods())>4) detail <- "This might take a while..."
      withProgress(message=msg, detail=detail, value=1, max=3, {
        hits$hits <- findSeedMatches( target(), selmods(),
                                      keepMatchSeq=input$keepmatchseq,
                                      minDist=input$minDist,
                                      shadow=input$shadow,
                                      BP=SerialParam(progressbar=TRUE) )
        if(length(hits$hits)>0){
          hits$hits$miRNA <- hits$hits$seed
          hits$hits$seed <- NULL
        }
      })
      if(input$subjet_type=="custom"){
        hits$target <- "custom sequence"
      }else{
        hits$target <- paste(input$gene, "-", seltx(),
                             ifelse(input$utr_only, "(UTR)",""))
      }
    })
    
    output$scan_target <- renderText({
      if(is.null(hits$target)) return(NULL)
      paste("Scan results in: ", hits$target)
    })
    
    output$hits_table <- renderDT({
      if(is.null(hits$hits)) return(NULL)
      h <- as.data.frame(hits$hits)
      h <- h[,setdiff(colnames(h), c("seqnames","width","strand") )]
      if(hits$nsel == 1) h$miRNA <- NULL
      h
    })
    
    output$manhattan <- renderPlot({
      if(is.null(hits$hits)) return(NULL)
      h <- as.data.frame(sort(hits$hits))
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

  }
}
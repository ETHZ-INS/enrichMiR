enrichMiR.server <- function(modlists){
  function(input, output, session){
    
    ##############################
    ## select collection
    
    updateSelectizeInput(session, "mirlist", choices=names(modlists))
    
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
    
    target <- reactive({ # target subject sequence
      if(input$subjet_type=="custom") return(input$customseq)
      # TO DO: fetch and return the selected transcript's sequence
      return(NULL)
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
      hits$hits <- NULL
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
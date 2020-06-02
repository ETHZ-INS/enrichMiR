enrichMiR.server <- function(modlists){
  function(input, output, session){
    updateSelectizeInput(session, "mirlist", choices=names(modlists))
    
    mods <- reactive({
      if(is.null(input$mirlist)) return(NULL)
      modlists[[input$mirlist]]
    })
    
    output$collection_summary <- renderPrint({
      if(is.null(mods)) return(NULL)
      #summary(mods)
      return(NULL)
    })
    
    observe({
      updateSelectizeInput(session, "mirnas", choices=names(mods()), server=TRUE)
      updateSelectizeInput(session, "mirna", choices=names(mods()), server=TRUE)
    })
    
    target <- reactive({
      if(input$subjet_type=="custom") return(input$customseq)
      # fetch and return the selected transcript's sequence
      return(NULL)
    })
    
    
    ###### begin miRNA-centric tab
    
    mod <- reactive({
      if(is.null(mods()) || is.null(input$mirna)) return(NULL)
      mods()[[input$mirna]]
    })
    
    output$modplot <- renderPlot({
      if(is.null(mod())) return(NULL)
      plotKdModel(mod())
    }, height=reactive(paste0(input$modplot_height,"px")))

  }
}
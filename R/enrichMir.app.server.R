#' enrichMiR.server
#'
#' @return A shiny server function
#' @export
#' @import ggplot2 DT GO.db
enrichMiR.server <- function(){
  library(DT)
  library(ggplot2)
  library(enrichMiR)
  library(GO.db)
  
  dtwrapper <- function(d, pageLength=25){
    datatable( d, filter="top", class="compact", extensions=c("Buttons","ColReorder"),
               options=list(pageLength=pageLength, dom = "fltBip", rownames=FALSE,
                            colReorder=TRUE, 
                            buttons=c('copy', 'csv', 'excel', 'csvHtml5') ) )
  }
 
  
# this first part doesn't work yet  
  
  
############################
#############################
  
  
  
   
  
  function(input, output, session){
    
    
    ##############################
    ## Select Species
    
    # this should be replaced by the TS object
    all_fams <- readRDS("/mnt/schratt/enrichMiR/data/20201006_Targetscan_Human_miRFamilies.rds")
    Species_fam <- reactive({
      # to be replaced with data()
      spec <- switch( input$species,
                      "Human" = 9606,
                      "Mouse" = 10090,
                      "Rat" = 10116,
                      "Custom - not yet" = 9606)
      fam <- all_fams[all_fams$`Species ID` == spec,]
      x <- fam$`Seed+m8`
      names(x) <- fam$`MiRBase ID`
      x
    })
    
    
    
    ##############################
    ## initialize expression info
    
    DEA <- reactive({ #initialize dea
      upFile <- input$dea_input
      #if (is.null(upFile)) return(NULL)
      if (is.null(upFile)){ ## TEMPORARY - BECAUSE I'M LAZY
        updf <- readRDS("/mnt/schratt/enrichMiR_benchmark/data/bartel_HEK.deaList.rds")[[1]]
      }else{
        updf <- read.csv(upFile$datapath, header = input$header, row.names=1)
      }
      enrichMiR:::.homogenizeDEA(updf)
    })
    
    
    Back <- reactive({ #initalize the background
      if(input$input_type == "dea"){
        row.names(DEA())
      }else{
        return(unlist(strsplit(input$background_genes, "\n")))
      }
    })
        
    miRNA <- reactive({ #initialize expressed miRNAS
      if (is.null(input$expressed_mirnas)) return(NULL)
      return(unlist(strsplit(input$expressed_mirnas, "\n")))
    })
        
    
    
    ##############################
    ## initialize reactive inputs
    
    # miRNA families
    observe(updateSelectizeInput(session, "mir_fam", choices=Species_fam(), server=TRUE))
    
    
    # Add GO Terms to input list
    
    GO_all <- as.data.frame(GO.db::GOTERM)
    GO_all_vec <- GO_all$go_id
    names(GO_all_vec) <- paste0(GO_all$go_id," (",GO_all$Term,")")
    GO_all_vec <- GO_all_vec[!duplicated(GO_all_vec)]
    updateSelectizeInput(session, "go_term", choices=GO_all_vec, server=TRUE)
    
    
    ##############################
    ## Initialize Genes of Interest
    
    Gene_Subset <- reactive({
      if(input$input_type == "dea"){
        if(is.null(input$dea_input)) return(NULL)
        d <- DEA()
        return(row.names(d)[d$FDR<input$binary_fdr])
      }else{
        if(input$GOI == "GOI_custom"){
          return(unlist(strsplit(input$genes_of_interest, "\n")))
        }else{
            Sp <- switch(input$species,
                                 "Human" = "Hs",
                                 "Mouse" = "Mm",
                                 "Rat" = "Rn",
                                 "Custom - not yet" = "Hs")
            ens <- switch(input$go_genes_format,
                                  "Ens" = TRUE,
                                  "GS" = FALSE)
            return(as.character(unlist(getGOgenes(go_ids = input$go_term, species = Sp,ensembl_ids = ens))))
          }
        }
    })
    
    output$GOI_nb <- renderText({
      genes <- Gene_Subset()
      if(is.null(genes)) return(NULL)
      paste(length(genes), " gene(s)")
    })
    
    ##############################
    ## Initialize target predictions
    
    TS <- reactive({
      # DO SOMETHING HERE
      # temporary, for testing:
      TS <- readRDS("/mnt/schratt/enrichMiR_benchmark/data/TargetScan_human.rds")
      colnames(TS)[1] <- "set"
      TS
    })
    
    
    ##############################
    ## CD plot
    
    output$cd_plot <- renderPlot({
      CDplot2(DEA(), TS(), setName=input$mir_fam, by = input$CD_type)
    })
    
    ##############################
    ## Enrichment analysis
    
    ER <- reactive({
      if(input$input_type == "dea"){
        if(is.null(DEA())) return(NULL)
        return(enrichMiR(DEA(), TS = TS(), tests=c("siteoverlap","areamir")))
      }else{
        if(is.null(Gene_Subset()) || is.null(Back())) return(NULL)
        return(testEnrichment(Gene_Subset(), TS(), background=Back(), tests=c("siteoverlap","areamir")))
      }
    })
    
    output$bubble_plot <- renderPlotly({
      if(is.null(ER())) return(NULL)
      ggplotly(enrichPlot(getResults(ER(), "areamir"), repel=FALSE))
    })    
    
  }
}





  
  
  



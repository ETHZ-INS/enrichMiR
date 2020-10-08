#' enrichMiR.server
#'
#' @return A shiny server function
#' @export
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
 
  
  
  
  
############################
#############################
  
  
  
   
  
enrichMiR.server <-function(input, output, session){
    
    
    ##############################
    ## Select Species
    
    Species_fam <- reactive({
      fam <- readRDS("data/20201006_Targetscan_Human_miRFamilies.rds")
      spec <- switch( input$species,
                      "Human" = 9606,
                      "Mouse" = 10090,
                      "Rat" = 10116,
                      "Custom - not yet" = 9606)
      fam[fam$`Species ID` == spec, c("miR family","Seed+m8","MiRBase ID")]
    })
    
    
    
    ##############################
    ## initialize expression info
    
    DEA <- reactive({ #initialize dea
      upFile <- input$dea_input
      if (is.null(upFile)) return(NULL)
      updf <- read.csv(upFile$datapath, header = input$header)
      updf
      })
    
    
    Back <- reactive({ #initalize the background
      if(input$upload_background=="upload"){
        return(unlist(strsplit(input$background_genes, "\n")))
      }else{
        DEA()[,1]
      }
    })
        
    miRNA <- reactive({ #initialize expressed miRNAS
      if (is.null(input$expressed_mirnas)) return(NULL)
      return(unlist(strsplit(input$expressed_mirnas, "\n")))
        })
        
    
    
    ##############################
    ## initialize reactive inputs
    
    # miRNA families
    # As long as I have a "data.frame" as input, one has to stay with 'selectInput'
    # >> Can be changed
    observe(updateSelectInput(session, "mir_fam", choices=Species_fam()))
    
    
    # Find GO Terms and update the input list
    GO <- eventReactive(input$find_go_button, {
      a <- suppressWarnings(enrichMiR::findGO(expr = as.character(input$find_go), ontology = as.character(input$ontology)))
      if(is.null(a)) return("Nothing Found")
      a
    })
    
    output$find_go_result <- renderPrint({
      if(is.null(GO())) return(NULL)
      head(GO(),20)
      })

    
    GO_all <- as.data.frame(GO.db::GOTERM)
    GO_all_vec <- paste0(GO_all$go_id," (",GO_all$Term,")")
    GO_all_vec <- GO_all_vec[!duplicated(GO_all_vec)]
    names(GO_all_vec) <- GO_all_vec

    GO_Input <- reactive({
      if(input$find_go_button == 0){
      return(GO_all_vec)
      }else{
      paste0(names(GO())," (",as.character(GO()),")")
      }
    })
    
    observe(updateSelectizeInput(session, "go_term", choices=GO_Input(),server = TRUE))
    
    
    
    ##############################
    ## Initialize Genes of Interest
    
    Gene_Subset <- eventReactive(input$enrich | input$colocalize, {#Get the genes of interest
      if(input$choose_gene_set == "dea_sets"){
        if(is.null(input$dea_input)) return(NULL)
        {d <- DEA()
        d[d[[3]] <= input$binary_fdr , 1]}
      }else{
          if(input$cust_enrichment_type == "custom_genes"){
            return(unlist(strsplit(input$genes_of_interest, "\n")))
          }else{
            a <- unlist(strsplit(input$go_term, " "))[1]
            Sp <- switch(input$species,
                                 "Human" = "Hs",
                                 "Mouse" = "Mm",
                                 "Rat" = "Rn",
                                 "Custom - not yet" = "Hs")
            ens <- switch(input$go_genes_format,
                                  "Ens" = TRUE,
                                  "GS" = FALSE)
            return(as.character(unlist(getGOgenes(go_ids = a,species = Sp,ensembl_ids = ens))))
          }
        }
    })
    
    
    observeEvent(input$enrich | input$colocalize,{
      print(head(Gene_Subset()))
    })
    
    ##############################
    ## Perform the enrichment analysis
    
    #if dea == cont, then cont test with dea
    #everywhen else use binary test with set = goi and back
    
    
    
    
    
  }  
  



shinyApp(enrichMiR.ui,enrichMiR.server) 
  


}





  
  
  



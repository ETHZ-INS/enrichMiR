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
    # MS>: the TS only contains repr. miRNAs, therefore I wanted this one to display all


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
    
    ## Include , and "" gsub
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



    ### here I have to subset for the expressed miRNAs
    
    #should this somehow be cached? Wanted to have it as download, because the files are so big
    TS_nonC <- reactive({
      sp <- switch(input$species,
                   "Human" = "human",
                   "Mouse" = "mouse",
                   "Rat" = "rat",
                   "Custom - not yet" = "human")
      TS_nonC <- .loadTargetscanSitesAll(species = sp)
      TS_nonC$sites <- TS_nonC$`Total num conserved sites` + TS_nonC$`Total num nonconserved sites`
      TS_nonC <- TS_nonC[,c("Gene Symbol","miRNA family","Cumulative weighted context++ score","sites")]
      colnames(TS_nonC) <- c("feature","set","score","sites")
      TS_nonC
    })
    
    TS_con <- reactive({
      TS_con <- switch(input$species,
                   "Human" = readRDS("/mnt/schratt/enrichMiR/data/20201006_Targetscan_Human_ConSites_Human.rds"),
                   "Mouse" = readRDS("/mnt/schratt/enrichMiR/data/20201006_Targetscan_Mouse_ConSites_Mouse.rds"),
                   "Rat" = readRDS("/mnt/schratt/enrichMiR/data/20201006_Targetscan_Mouse_ConSites_Rat.rds"),
                   "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR/data/20201006_Targetscan_Human_ConSites_Human.rds"))
      TS_con$sites <- TS_con$`Total num conserved sites`
      TS_con <- TS_con[,c("Gene Symbol","miRNA family","Cumulative weighted context++ score","sites")]
      colnames(TS_nonC) <- c("feature","set","score","sites")
      TS_con
    })
    
    Scan_Mir <- reactive({
      detail <- "not yet..."
    })
    
    RBPs <- reactive({
      detail <- "not yet..."
    })
    
    #mirTarbas
    
    Custom <- reactive({
      detail <- "not yet..."
    })
      
    EN_Object <- reactive({
      EN_Object <- switch(input$collection,
                          "scanMir miRNA BS" = Scan_Mir(),
                          "Targetscan conserved miRNA BS" = TS_con(),
                          "Targetscan all miRNA BS"= TS_nonC(),
                          "CISBP RBP motif sites" = RBPs(),
                          "Custom - not yet" = Custom())
      
      if (is.null(input$expressed_mirnas)) return(EN_Object)
      if(input$collection == "Targetscan conserved miRNA BS" | input$collection == "Targetscan all miRNA BS"){
        EN_Object <- EN_Object[EN_Object$set %in% 
                                 
                                 
      #############
      #continue here with the miRNA filtering
      #built the miRNA family + name + motif data.frame up where I initalize miRNAs 
                                 
                                 
                       
      })



    
    ##############################
    ## CD plot
    
    output$cd_plot <- renderPlot({
      CDplot2(DEA(), EN_Object(), setName=input$mir_fam, by = input$CD_type)
    })
    
    ##############################
    ## Enrichment analysis
    
    ER <- reactive({
      if(input$collection == "Targetscan all miRNA BS") detail <- "This will take a while..."
      if(input$input_type == "dea"){
        if(is.null(DEA())) return(NULL)
        return(enrichMiR(DEA(), TS = EN_Object(), tests=c("siteoverlap","areamir")))
      }else{
        if(is.null(Gene_Subset()) || is.null(Back())) return(NULL)
        return(testEnrichment(Gene_Subset(), EN_Object(), background=Back(), tests=c("siteoverlap","areamir")))
      }
    })
    
    output$bubble_plot <- renderPlotly({
      if(is.null(ER())) return(NULL)
      ggplotly(enrichPlot(getResults(ER(), "areamir"), repel=FALSE))
    })    
    
  }
}





  
  
  



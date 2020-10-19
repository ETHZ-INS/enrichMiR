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
        
    miRNA_exp <- reactive({ #initialize expressed miRNAS
      miRNA_up <- input$miRNA_exp_input
      if (is.null(upFile)){return(NULL)
      }else{
        mirup <- read.csv(miRNA_up$datapath, header = input$header_mir, row.names=1)
      }
      mirup <- mirup[order(mirup[[2]]),]
      mirup <- mirup[1:((input$mir_cut_off/100)*nrow(mirup)),]
      mirup
    })
        
    
    miRNA_exp_list <-reactive({ #initialize custom list of expressed miRNAS 
      if (is.null(input$expressed_mirnas)) return(NULL)
      return(unlist(strsplit(input$expressed_mirnas, "\n")))
    })
    
    ##############################
    ## initialize reactive inputs
  
    # Add GO Terms to input list
    GO_all <- as.data.frame(GO.db::GOTERM)
    GO_all_vec <- GO_all$go_id
    names(GO_all_vec) <- paste0(GO_all$go_id," (",GO_all$Term,")")
    GO_all_vec <- GO_all_vec[!duplicated(GO_all_vec)]
    updateSelectizeInput(session, "go_term", choices=GO_all_vec, server=TRUE)
    
    output$GOI_nb <- renderText({
      genes <- Gene_Subset()
      if(is.null(genes)) return(NULL)
      paste(length(genes), " gene(s)")
    })
    
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
    

    
    
    ##############################
    ## Load En_Objects based on species
    
    
    ### Include here the miRNA filtering
    
 #    if (is.null(input$exp_mirna_list)){
 #      if(is.null(input$exp_mirna_file)){
 #        return(EN_Object)
 #      }else{
 #        #include expression to family metadata and filter for exp. miRNAs  
 #      }
 #    }else{
 #      #filter for expressed miRNAs in list
 #    }
 # })
    
    TS_nonC <- reactive({
      if(input$genes_format == "Ens"){
        TS_nonC <- switch(input$species,
                          "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Human_AllSites_human_InE_Fam.rds"),
                          "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Mouse_AllSites_mouse_InE_Fam.rds"),
                          "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Mouse_AllSites_rat_InE_Fam.rds"),
                          "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Human_AllSites_human_InE_Fam.rds"))
      }
      if(input$genes_format == "GS"){
        TS_nonC <- switch(input$species,
                          "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Human_AllSites_human_InS_Fam.rds"),
                          "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Mouse_AllSites_mouse_InS_Fam.rds"),
                          "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Mouse_AllSites_rat_InS_Fam.rds"),
                          "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Human_AllSites_human_InS_Fam.rds"))
      }
      TS_nonC
    })
    
    TS_con <- reactive({
      if(input$genes_format == "Ens"){
        TS_con <- switch(input$species,
                         "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Human_ConSites_human_InE_Fam.rds"),
                         "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Mouse_ConSites_mouse_InE_Fam.rds"),
                         "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Mouse_ConSites_rat_InE_Fam.rds"),
                         "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Human_ConSites_human_InE_Fam.rds"))
      }
      if(input$genes_format == "GS"){
        TS_con <- switch(input$species,
                         "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Human_ConSites_human_InS_Fam.rds"),
                         "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Mouse_ConSites_mouse_InS_Fam.rds"),
                         "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Mouse_ConSites_rat_InS_Fam.rds"),
                         "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Human_ConSites_human_InS_Fam.rds"))
      }
      print(head(TS_con))
      TS_con
    })
    
    Scan_Mir <- reactive({
      detail <- "not yet..."
    })
    
    RBPs <- reactive({
      detail <- "not yet..."
    })
    
    miRTarBase <- reactive({
      detail <- "not yet..."
    })
    
    Custom <- reactive({
      detail <- "not yet..."
    })
    
    
    
    
    
    
    ##############################
    ## Initialize target predictions

    
    EN_Object <- reactive({
      switch(input$collection,
             "scanMir miRNA BS"=switch(input$species,
                                       ))
      # switch(input$collection,
      #         "scanMir miRNA BS" = Scan_Mir(),
      #         "Targetscan conserved miRNA BS" = TS_con(),
      #         "Targetscan all miRNA BS"= TS_nonC(),
      #         "CISBP RBP motif sites" = RBPs(),
      #         "miRTarBase" = miRTarBase(),
      #         "Custom - not yet" = Custom())
      readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201017_Targetscan_Mouse_ConSites_rat_InE_Fam.rds")
    })
                                 
   
    
   
    ##############################
    ## CD plot
      
    ### Include here the S4 metadata!!
    observe({
      print(head(EN_Object()))
      if(!is.null(EN_Object())){
        message("EN not null")
        if(!is.null(m <- metadata(EN_Object())$families)){
          ids <- m[["Seed+m8"]]
          names(ids) <- [["MiRBase ID"]]
          updateSelectizeInput(session, "mir_fam", choices=unique(ids), server=TRUE)
        }
      }
    }) 
              
        
              
    
    output$cd_plot <- renderPlot({
      if(is.null(input$mir_fam) || input$mir_fam=="") return(NULL)
      CDplot2(DEA(), EN_Object(), setName=input$mir_fam, by = input$CD_type)
    })
    
    ##############################
    ## Enrichment analysis
    
    ## include options from continuous tab??
    ## Decide between continuous and binary test already in ER or only in plotting / data.table
    
    
    ER <- reactive({
      if(input$collection == "Targetscan all miRNA BS") detail <- "This will take a while..."
      if(input$input_type == "dea"){
        if(is.null(DEA())) return(NULL)
        return(enrichMiR(DEA(), TS = EN_Object(), tests=c("siteoverlap","areamir"), minSize=input$minsize))
      }else{
        if(is.null(Gene_Subset()) || is.null(Back())) return(NULL)
        return(testEnrichment(Gene_Subset(), EN_Object(), background=Back(), tests=c("siteoverlap","areamir"),minSize=input$minsize))
      }
    })
    
    output$bubble_plot <- renderPlotly({
      if(is.null(ER())) return(NULL)
      
      
      ggplotly(enrichPlot(getResults(ER(), "areamir"), repel=FALSE))
    })    
    
  }
}





  
  
  



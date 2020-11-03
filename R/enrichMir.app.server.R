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
  
  dtwrapper <- function(d, pageLength=25, hide_cols){
    datatable( d, filter="top", class="compact", extensions=c("Buttons","ColReorder"),
               options=list(pageLength=pageLength, dom = "fltBip", rownames=FALSE,
                            colReorder=TRUE,
                            columnDefs = list(list(visible=FALSE, 
                                                   targets=na.omit(match(hide_cols, colnames(d))))),
                            buttons=c('copy', 'csv', 'excel', 'csvHtml5') ) )
  }
  
  trimInputList <- function(x){
    x <- unlist(strsplit(gsub(",|;|\\t|\\r|\\s","\n",x),"\n"))
    x <- unique(x[x!=""])
    if(length(x)==0) return(NULL)
    x
  }
  
  

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
        if(is.null(DEA())) return(NULL)
        row.names(DEA())
      }else{
        return(trimInputList(input$background_genes))
      }
    })
        
    miRNA_exp <- reactive({
      if(isTRUE(getOption("shiny.testmode"))) print("miRNA_exp")
      if(is.null(input$expressed_mirna_box)) return(NULL)
      if(input$expressed_mirna_box=="Custom Set"){
        if(is.null(input$expressed_mirnas) || input$expressed_mirnas=="") return(NULL)
        return(trimInputList(input$expressed_mirnas))
      }
      if(is.null(input$miRNA_exp_input)) return(NULL)
      mirup <- read.csv(input$miRNA_exp_input$datapath, header = input$header_mir, row.names=1)
      mirup <- mirup[order(mirup[[2]]),]
      colnames(mirup)[1] <- "name"
      colnames(mirup)[2] <- "expression"
      mirup <- mirup[1:((input$mir_cut_off/100)*nrow(mirup)),]
      return(mirup)
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
      if(isTRUE(getOption("shiny.testmode"))) print("Gene_Subset")
      if(is.null(input$input_type)) return(NULL)
      if(input$input_type == "dea"){
        if(is.null(input$dea_input)) return(NULL)
        d <- DEA()
        return(row.names(d)[d$FDR<input$binary_fdr])
      }else{
        if(input$GOI == "GOI_custom"){
          return(trimInputList(input$genes_of_interest))
        }else{
            Sp <- switch(input$species,
                                 "Human" = "Hs",
                                 "Mouse" = "Mm",
                                 "Rat" = "Rn",
                                 "Custom - not yet" = "Hs")
            ens <- switch(input$genes_format,
                                  "Ens" = TRUE,
                                  "GS" = FALSE)
            return(as.character(unlist(getGOgenes(go_ids = input$go_term, species = Sp,ensembl_ids = ens))))
          }
        }
    })
    

  
    ##############################
    ## Initialize target predictions

    
    EN_Object <- reactive({
      if(is.null(EN_Object)) return(NULL)
      if(isTRUE(getOption("shiny.testmode"))) print("EN_Object")
      switch(input$collection,
            #Loading everywhere the conserved files for now, except for the "all_Sites" object
             "scanMir miRNA BS" = switch(input$species,
                                                   "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds"),
                                                   "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_mouse.rds"),
                                                   "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_rat.rds"),
                                                   "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds")),
             "Targetscan conserved miRNA BS" = switch(input$species,
                                                        "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds"),
                                                        "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_mouse.rds"),
                                                        "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_rat.rds"),
                                                        "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds")),
             "Targetscan all miRNA BS" = switch(input$species,
                                                       "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_AllSites_human.rds"),
                                                       "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_AllSites_mouse.rds"),
                                                       "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_AllSites_rat.rds"),
                                                       "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_AllSites_human.rds")),
             "CISBP RBP motif sites" = switch(input$species,
                                                         "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds"),
                                                         "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_mouse.rds"),
                                                         "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_rat.rds"),
                                                         "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds")),
             "miRTarBase" = switch(input$species,
                                              "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds"),
                                              "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_mouse.rds"),
                                              "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_rat.rds"),
                                              "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds")),
             "Custom - not yet" =  switch(input$species,
                                          "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds"),
                                          "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_mouse.rds"),
                                          "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Mouse_ConSites_rat.rds"),
                                          "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201020_Targetscan_Human_ConSites_human.rds")))
    })

                                 
   
    
   
    ##############################
    ## CD plot
    
    observe({
      if(!is.null(EN_Object())){
        if(!is.null(m <- metadata(EN_Object())$families)){
          updateSelectizeInput(session, "mir_fam", choices=m, server=TRUE)
        }
      }
    }) 
              
    output$cd_plot <- renderPlot({
      if(is.null(input$mir_fam) || input$mir_fam=="") return(NULL)
      CDplot2(DEA(), EN_Object(), setName=input$mir_fam, by = input$CD_type)
    })
    
    ##############################
    ## Enrichment analysis
    
    ViewTest <- reactive({
      switch(input$test_type,
             "binary" = switch(input$up_down,
                               ".up" = "siteoverlap.up",
                               ".down" = "siteoverlap.down"),
             "continous" = "areamir",
             "advanced" = input$view_test)
    })
    
    output$test_info <- renderPrint({ # print the test
      if(is.null(ViewTest()) || ViewTest()=="") return("")
      out <- capture.output(ViewTest())
      cat(out)
    })
      
      
    # We would like the the overlap names in the result table
    
    # We would like to have the mir_names in the result table
    
    ER <- eventReactive(input$enrich, {
      if(is.null(input$input_type) || is.null(EN_Object())) return(NULL)
      if(isTRUE(getOption("shiny.testmode"))) print("ER")
      mirexp <- miRNA_exp()
      if(!is.null(mirexp) && is.null(dim(mirexp)))
        mirexp <- data.frame(row.names=mirexp, miRNA=mirexp)
      
      standard_tests <- c("siteoverlap","areamir")
      tests <- c(standard_tests, input$tests2run)
      msg <- paste0("Performing enrichment analyses with the following tests: ",paste(tests,collapse=", "))
      message(msg)
      if(input$input_type == "dea"){
        if(is.null(DEA())) return(NULL)
        sig <- DEA()
        bg <- NULL
      }else{
        if(is.null(Gene_Subset()) || is.null(Back())) return(NULL)
        sig <- Gene_Subset()
        bg <- Back()
      }
      if(length(sig)==0) return(NULL)
      detail <- NULL
      if(input$collection == "Targetscan all miRNA BS") detail <- "This will take a while..."
      withProgress(message=msg, detail=detail, value=1, max=3, {
      testEnrichment(sig, EN_Object(), background=bg, sets.properties = mirexp,
                     tests=tests, minSize=input$minsize, th.FDR=input$dea_sig_th)
      })
    })
    
    observe({
      if(!is.null(ER())) updateSelectInput(session, "view_test", choices=c("",names(ER())))
    })

    
    output$bubble_plot <- renderPlotly({
      if(is.null(ER())) return(NULL)
      test <- ViewTest()
      if(is.null(ViewTest()) || ViewTest()==""){
        er <- getResults(ER(), getFeatures=FALSE, flatten=TRUE)
        er$FDR <- er$FDR.geomean
        er$enrichment <- rowMeans(er[,grep("nrichment|beta|coefficient",colnames(er)),drop=FALSE],na.rm=TRUE)
      }else{
        er <- getResults(ER(), test=test, getFeatures=FALSE, flatten=TRUE)
      }
      ggplotly(enrichPlot(er, repel=FALSE, label.sig.thres = input$label.sig.thres,
                          min.enr.thres = input$label.enr.thres, maxLabels = input$label_n ))
    })
    
    # columns2hide <- reactive({
    #   show_standard <- c("enrichment","pvalue","FDR")
    #   if(is.null(input$columns2show)) return(setdiff(colnames(output$hits_table),show_standard))
    #   show_add <- input$columns2show
    #   if(grepl("genes",show_add)){
    #     show_add["genes"] <- NULL
    #     show_add <- c(show_add,"genes.down","genes.up")
    #   }
    #   show <- c(show_standard,show_add)
    #   return(setdiff(colnames(output$hits_table),show))
    # })
    
    output$hits_table <- renderDT({ # prints the current hits
      if(is.null(ER())) return(NULL)
      test <- ViewTest()
      if(is.null(ViewTest()) || ViewTest()=="") test <- NULL
      rr <- getResults(ER(), test=test, flatten=TRUE)
      show_standard <- c("enrichment","pvalue","FDR")
      if(is.null(input$columns2show)){
        columns2hide <- setdiff(colnames(rr),show_standard)
      }else{
        show_add <- input$columns2show
          if(any(show_add %in% "genes")){
            show_add <- c(show_add,"genes.down","genes.up")
            show_add <- show_add[!show_add %in% "genes"]
            show <- c(show_standard,show_add)
          }else{
            show <- c(show_standard,show_add)
          }
        columns2hide <- setdiff(colnames(rr),show)
      }
      return(dtwrapper(rr,hide_cols = columns2hide))
    })
    
    output$dl_hits <- downloadHandler(
      filename = function() {
        if(is.null(ER())) return(NULL)
        fn <- paste0("EnrichMir_hits_",ViewTest(),"_",Sys.Date(),".csv")
        fn
      },
      content = function(con) {
        if(is.null(ER())) return(NULL)
        test <- ViewTest()
        h <- getResults(ER(), test=test, flatten=TRUE)
        write.csv(h, con)
      }
    )
      
    
    if(isTRUE(getOption("shiny.testmode"))) print("END LOADING")
    
  }
}





  
  
  



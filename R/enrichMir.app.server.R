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
        if(input$background_genes == "") return(NULL)
        return(unlist(strsplit(input$background_genes, "\n")))
      }
    })
        
    miRNA_exp_file <- reactive({ #initialize expressed miRNAS
      miRNA_up <- input$miRNA_exp_input
      if (is.null(upFile)){return(NULL)
      }else{
        mirup <- read.csv(miRNA_up$datapath, header = input$header_mir, row.names=1)
      }
      mirup <- mirup[order(mirup[[2]]),]
      colnames(mirup)[1] <- "name"
      colnames(mirup)[2] <- "expression"
      mirup <- mirup[1:((input$mir_cut_off/100)*nrow(mirup)),]
      mirup
    })
        
    
    miRNA_exp_list <-reactive({ #initialize custom list of expressed miRNAS 
      if (input$expressed_mirnas =="") return(NULL)
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
          if(input$genes_of_interest == "") return(NULL)
          return(unlist(strsplit(input$genes_of_interest, "\n")))
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
      eno <- switch(input$collection,
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
      if(is.null(eno)) return(NULL)
      
      #This part here is for the miRNA filterin and creating family row.names in the metadata
      if(!is.null(metadata(eno)$families)){
      m <- metadata(eno)$families
      if(is.null(miRNA_exp_list())){
        if(is.null(miRNA_exp_file())){
          ll <- CharacterList(lapply(split(names(m), m),unique))
          print(head(ll))
          ll <- vapply(ll,FUN.VALUE = character(length = 1),FUN=function(x){ paste(x, collapse = ", ") })
          df <- as.data.frame(ll,row.names = names(ll))
          metadata(eno)$families <- df
          return(eno)
        }else{
          # filter here for the uploaded miRNAs and assign the expression value (at the moment as sum of the family)
          m <- as.data.frame(fam = m, name = names(m))
          m <- merge(m,miRNA_exp(),by = "name")
          ll <- CharacterList(lapply(split(m$name, m$fam),unique))
          ll <- vapply(ll,FUN.VALUE = character(length = 1),FUN=function(x){ paste(x, collapse = ", ") })
          df <- as.data.frame(ll,row.names = names(ll))
          df$fam <- row.names(df)
          m <- m %>% group_by(fam = m$fam) %>%
            summarize(exp = sum(expression))
          df <- merge(df,m,by = "fam" )
          metadata(eno)$families <- df
          return(eno)
        }
      }else{
        #filter here for the miRNAs that are pasted
        m <- m[names(m) %in% input$exp_mirna_list]
        m <- droplevels(m)
        ll <- CharacterList(lapply(split(names(m), m),unique))
        ll <- vapply(ll,FUN.VALUE = character(length = 1),FUN=function(x){ paste(x, collapse = ", ") })
        df <- as.data.frame(ll,row.names = names(ll))
        metadata(eno)$families <- df
        return(eno)
      }

      }


      })

                                 
   
    
   
    ##############################
    ## CD plot
    
    # Doesn't work yet. Error seems to lie above in the metadata assignment
    
      
    observe({
      if(!is.null(EN_Object())){
        if(!is.null(m <- metadata(EN_Object())$families)){
          ids <- names(m)
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
    
    # I would do the Ensembl / Gene Symbol filtering in here since we need to do it seperate for DEAs and Subset / Back
    # >> I tried but this doesn't work yet, doesn't convert the feature names to Ensembl
    
    # Include the observeEvent
    
    # We would like the the overlap names in the result table
    
    # We would like to have the mir_names in the result table
    
    
    ER <- reactive({
      if(input$collection == "Targetscan all miRNA BS") detail <- "This will take a while..."
      if(input$input_type == "dea"){
        if(is.null(DEA())) return(NULL)
          # if(!is.null(metadata(EN_Object())$feature.synonyms) &&
          #  length(w <- which(row.names(DEA()) %in% names(metadata(EN_Object())$feature.synonyms))>0)){
          #    g <- row.names(DEA())
          #    g[w] <- metadata(EN_Object())$feature.synonyms[g[w]]
          #    g <- g[w <- !duplicated(g)]
          #    DEA() <- DEA()[w,]
          #    row.names(DEA()) <- g
          # }
        print(head(EN_Object()))
        return(enrichMiR(DEA(), TS = EN_Object(), tests=c("siteoverlap","areamir"), minSize=input$minsize,th.FDR = input$dea_sig_th))
      }else{
        if(is.null(Gene_Subset()) || is.null(Back())) return(NULL)
        # if(input$genes_format == "ENS"){
        #   g <- Back()
        #   g[w] <- metadata(TS)$feature.synonyms[g[w]]
        #   g <- g[w <- !duplicated(g)]
        #   
        return(testEnrichment(Gene_Subset(), EN_Object(), background=Back(), tests=c("siteoverlap"),minSize=input$minsize))
      }
    })
    
    #include width or heigts, whichever was working with renderplotly
    
    output$bubble_plot <- renderPlotly({
      if(is.null(ER())) return(NULL)
      if(input$test_type == "binary"){
        if(input$input_type == "dea"){
          a <- paste0("siteoverlap",input$up_down)
        }else{
          a < "siteoverlap"
        }
        b <- getResults(ER(), a)
        sf <- "overlap"
      }else{
        a <- "areamir"
        b <- getResults(ER(), a)
        sf <- "fam_size"
      }
      ggplotly(enrichPlot(b, size.field = sf, repel=FALSE, label.sig.thres = input$label.sig.thres,
                          min.enr.thres = input$label.enr.thres,maxLabels = input$label_n ))
    })    
    
    output$hits_table <- renderDT({ # prints the current hits
      if(is.null(ER)) return(NULL)
      if(input$test_type == "binary"){
        if(input$input_type == "dea"){
          a <- paste0("siteoverlap",input$up_down)
        }else{
          a < "siteoverlap"
        }
      }else{
        a <- "areamir"
      }
      h <- as.data.frame(getResults(ER(),test = a))
      h$miRNA.family <- row.names(h)
      dtwrapper(h)
    })
    
    
    
  }
}





  
  
  



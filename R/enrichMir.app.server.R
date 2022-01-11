#' enrichMiR.server
#' 
#' @param bData A named, nested list, the first level being species, the second
#' different binding annotations (see vignette for format details)
#' @param logCallsFile Optional path to a file when the date of enrichment calls
#' will be logged.
#'
#' @return A shiny server function
#' @export
#' @import ggplot2 DT GO.db
#' @importFrom rintrojs introjs
#' @importFrom shinyjqui jqui_resizable
#' @importFrom shinyjs hideElement showElement
#' @importFrom plotly renderPlotly ggplotly event_data
enrichMiR.server <- function(bData=NULL, logCallsFile=NULL){
  library(DT)
  library(ggplot2)
  library(enrichMiR)
  library(GO.db)
  library(rintrojs)
  library(shinyjs)
  library(shinyjqui)
  
  baseDataPath="/mnt/schratt/enrichMiR_data/"
  if(is.null(bData)) bData <- lapply(list(
    Human=c(
      "Targetscan conserved miRNA BS"=
        "Targetscan/20211124_Targetscan8_Human_ConPred_human.rds",
      "Targetscan all miRNA BS"=
        "Targetscan/20211124_Targetscan8_Human_AllPred_human.rds",
      "miRTarBase"=
        "miRTarBase/miRTarBase8.human.rds",
      "scanMiR (GRCh38)"="scanMiR/scanMiR_GRCh38_merged.rds",
      "oRNAment"=
        "oRNAment/human_coding_MSS05_gene_3UTR_pred.rds"),
    Mouse=c(
      "Targetscan conserved miRNA BS"=
        "Targetscan/20211124_Targetscan8_Mouse_ConPred_mouse.rds",
      "Targetscan all miRNA BS"=
        "Targetscan/20211124_Targetscan8_Mouse_AllPred_mouse.rds",
      "miRTarBase"=
        "miRTarBase/miRTarBase8.mouse.rds",
      "scanMiR (GRCm38)"="scanMiR/scanMiR_GRCm38_merged.rds",
      "oRNAment"="oRNAment/mouse_coding_MSS05_gene_3UTR_pred.rds"),
    Rat=c(
      "Targetscan conserved miRNA BS (from mouse)"=
        "Targetscan/20211124_Targetscan8_Mouse_ConPred_rat.rds",
      "Targetscan all miRNA BS (from mouse)"=
        "Targetscan/20211124_Targetscan8_Mouse_AllPred_rat.rds",
      "miRTarBase"=
        "miRTarBase/miRTarBase8.rat.rds",
      "scanMiR (Rnor6)"="scanMiR/scanMiR_rnor6_merged.rds")
    ), function(x) setNames(paste0(baseDataPath,x), names(x)))
  
  dtwrapper <- function(d, pageLength=25, hide_cols){
    datatable( d, filter="top", class="compact", extensions=c("Buttons","ColReorder"),
               options=list(pageLength=pageLength, dom = "fltBip", rownames=FALSE,
                            colReorder=TRUE,
                            columnDefs=list(
                              list(visible=FALSE, 
                                   targets=na.omit(match(hide_cols, colnames(d))))),
                            buttons=c('copy', 'csv', 'excel') ) )
  }
  
  trimInputList <- function(x){
    x <- unlist(strsplit(gsub(",|;|\\t|\\r|\\s","\n",x),"\n"))
    x <- unique(x[x!=""])
    x <- x[!is.na(x)]
    if(length(x)==0) return(NULL)
    x
  }
  
  setTheme <- function(p, theme=NULL){
    if(is.null(theme)) return(p)
    tryCatch({
      if(theme %in% c("classic2", "pubr", "pubclean")){
        theme <- getFromNamespace(paste0("theme_",theme), "ggpubr")
      }else{
        theme <- getFromNamespace(paste0("theme_",theme), "ggplot2")
      }
      p + theme()
    }, error=function(e){ warning(e); p })
  }

 function(input, output, session){

   bDataLoaded <- lapply(bData, as.list)
   updateSelectInput(session, "species", choices=names(bData))
   
   ##############################
   ## Introduction & help
   ## See app.intro.R for the actual content
   
   startIntro <- function(session){
     introjs(session, options=list(steps=.getAppIntro(), "nextLabel"="Next",
                                   "prevLabel"="Previous"),
             events=list(onbeforechange=readCallback("switchTabs")))
   }
   
   observeEvent(input$helpLink, startIntro(session))
   observeEvent(input$helpBtn, startIntro(session))

   observeEvent(input$overview, showModal(.getHelpModal("overview")))
   observeEvent(input$help_collections, showModal(.getHelpModal("collections")))
   observeEvent(input$help_collections2, showModal(.getHelpModal("collections")))
   observeEvent(input$help_enrichplot, showModal(.getHelpModal("enrichplot")))
   observeEvent(input$help_cdplot, showModal(.getHelpModal("cdplot")))
   observeEvent(input$help_cdplot2, showModal(.getHelpModal("cdplot")))
   observeEvent(input$help_tests, showModal(.getHelpModal("tests")))
   observeEvent(input$help_tests2, showModal(.getHelpModal("tests2")))
   observeEvent(input$brow_comp, showModal(.getHelpModal("browsercompatibility")))
   observeEvent(input$help_testsadvanced, showModal(.getHelpModal("testsadvanced")))
   observeEvent(input$help_deaformat, showModal(.getHelpModal("deaformat")))
   
   output$testsummary <- renderTable(.getTestsTable())
   output$browsercomp <- renderTable(.getBrowserCompTable())
   
   runjs("
   $('.box').on('click', '.box-header h3', function() {
     $(this).closest('.box')
     .find('[data-widget=collapse]')
     .click();
   });
   ")
   
   ##############################
   ## Menu items
   
   output$menu_species <- renderMenu({
     menuItem( "Species and miRNAs", tabName="tab_species", 
               expandedName="menu_species", icon=icon("folder-open"), 
               badgeLabel=input$species, badgeColor="light-blue")
   })
   output$menu_input <- renderMenu({
     if(input$input_type=="dea" && !is.null(DEA()))
       return(menuItem("Input genes/DEA", tabName="tab_input", badgeLabel="DEA",
                       badgeColor="light-blue", icon=icon("file-alt"),
                       expandedName="menu_input"))
     if(input$input_type!="dea" && length(Gene_Subset())>1 && 
        length(Back())>1)
       return(menuItem("Input genes/DEA", tabName = "tab_input", 
                       badgeLabel=paste0("set(",length(Gene_Subset()),
                                         "/",length(Back()),")"),
                       badgeColor="light-blue", icon=icon("file-alt"),
                       expandedName="menu_input"))
     menuItem("Input genes/DEA", tabName = "tab_input", badgeLabel="none", 
              badgeColor="red", icon=icon("file-alt"), expandedName="menu_input")
   })
   output$menu_enrich <- renderMenu({
     if((input$input_type=="dea" && !is.null(DEA())) ||
        (input$input_type!="dea" && length(Gene_Subset())>1 && length(Back())>1))
       return(menuItem("Enrichment analysis", tabName = "tab_enrich", 
                       icon=icon("rocket"), expandedName="menu_enrich"))
     menuItem("Enrichment analysis", tabName="tab_enrich", badgeLabel="(no input!)",
              badgeColor="red", icon=icon("times-circle"), expandedName="menu_enrich")
   })
   output$menu_cdplot <- renderMenu({
     if(is.null(DEA()) || input$input_type!="dea")
       return(menuItem("CD Plot", tabName="tab_cdplot", badgeLabel="(requires DEA)",
                       badgeColor="red", icon=icon("times-circle"),
                       expandedName="menu_cdplot"))
     menuItem("CD Plot", tabName="tab_cdplot", icon=icon("poll"),
              expandedName="menu_cdplot")
   })
   
   observe({
     if(!is.null(input$species))
       updateSelectInput(session, "collection", choices=names(bData[[input$species]]))
   })
   observe({
     if(!is.null(input$collection)){
       if(grepl("oRNAment",input$collection)){
         hideElement("exprMirs_box")
         hideElement("columns2show")
       }else{
         showElement("exprMirs_box")
         showElement("columns2show")
       }
       hideElement("resultsbox")
     }
   })
   
    ##############################
    ## initialize expression info
  
    DEA <- reactiveVal(NULL)
   
    observeEvent(input$dea_input, {
      upFile <- input$dea_input
      if (is.null(upFile)) return(NULL)
      updf <- tryCatch({
        .homogenizeDEA(data.table::fread(upFile$datapath))
      }, error=function(e){
        showModal(modalDialog(
          title = "There was an error reading your file:",
          tags$pre(as.character(e)),
          tags$p("Be sure to format your DEA table correctly:"),
          tags$p("Provide ENSEMBL_ID or Gene Symbol as identifier in the ",
                 tags$b("first")," column, as well as logFC-values and 
                 FDR-values. The output formats of most RNAseq DEA packages 
                 should be automatically recognized.")
        ))
        NULL
      })
      hideElement("resultsbox")
      if(!is.null(updf) && nrow(updf)>1) DEA(updf)
    })
    observeEvent(input$example_dea, {
      data(exampleDEA, package="enrichMiR")
      if(input$species != "Human"){
        row.names(exampleDEA) <-
          tools::toTitleCase(tolower(row.names(exampleDEA)))
      }
      hideElement("resultsbox")
      DEA(exampleDEA)
    })
    
    output$dea_res <- renderUI({
      if(is.null(DEA()))
        return(tags$span(icon("exclamation-triangle"), "No valid file uploaded",
                         style="font-weight: bold; font-color: red;"))
      tagList( tags$span(icon("check-circle"), "Valid file",
                         style="font-weight: bold; font-color: forestgreen;"),
               tags$p(nrow(DEA())," features, of which ", 
                      sum(DEA()$FDR<=input$dea_sig_th), " have a significance
                      below or equal to the selected FDR threshold (",
                      input$dea_sig_th, ")."))
    })
    
    output$enrichbtn <- renderUI({
      if( (input$input_type == "dea" && is.null(DEA())) ||
          (input$input_type != "dea" && (
            is.null(Gene_Subset()) || length(Gene_Subset())<2 ||
            is.null(Back()))) )
        return(tags$span(icon("exclamation-triangle"), "No valid gene/DEA input!",
                         style="font-weight:bold; font-color:red; font-size:120%;"))
      actionButton(inputId="enrich", "Run enrichMir!", icon = icon("search"))
    })
    
    ## Include , and "" gsub
    Back <- reactive({ #initalize the background
      hideElement("resultsbox")
      if(input$input_type == "dea"){
        if(is.null(DEA())) return(NULL)
        return(row.names(DEA()))
      }else{
        b <- trimInputList(input$background_genes)
        return(gsub("\\..*","",b))
      }
    })
    
    output$mirexp_preset <- renderUI({
      if(tolower(input$species) == "human"){
        return(tagList(
          tags$h4("From the human microRNAome package"),
          selectizeInput("mirexp_human", "Select tissue/celltype:", 
                         choices=getHumanMirExp()),
          tags$p("Source: ", 
                 tags$a(href="http://genome.cshlp.org/content/27/10/1769",
                        "McCall et al., NAR 2017")),
          
        ))
      }else if(tolower(input$species) == "mouse"){
        return(tagList(
          tags$h4("From published mouse profiles"),
          selectizeInput("mirexp_mouse", "Select tissue/celltype:",
                         choices=getMouseMirExp()),
          tags$p("Source: ",
                 tags$a(href="https://doi.org/10.1093/nar/gkaa323",
                        "Kern et al., NAR 2020"), " and ",
                 tags$a(href="https://doi.org/10.1016/j.neuron.2011.11.010",
                        "He et al., Neuron 2012")),
          tags$p("Expression is given in logCPM, i.e. log-transformed ",
                 "counts-per-million reads"),
          tags$p("Note that this quantification is not hairpin-specific,",
                 " but at the precursor level.")))
          
      }
      tags$p("Preset miRNA expression profiles are only available for human and mouse.")
    })
        
    miRNA_exp <- reactive({
      if(isTRUE(getOption("shiny.testmode"))) print("miRNA_exp")
      if(is.null(input$expressed_mirna_box)) return(NULL)
      if(input$expressed_mirna_box=="Custom Set"){
        if(is.null(input$exp_mirna_list) || input$exp_mirna_list=="") return(NULL)
        return(trimInputList(input$exp_mirna_list))
      }
      if(input$expressed_mirna_box=="Use preset expression profile"){
        if(tolower(input$species) == "human"){
          if(is.null(input$mirexp_human) || 
             is.null(x <- getHumanMirExp(input$mirexp_human))) return(NULL)
          x <- head(sort(x,decreasing=TRUE),round((input$mir_cut_off2/100)*length(x)))
          return(data.frame(row.names=names(x), expression=as.numeric(x)))     
        }else if(tolower(input$species) == "mouse"){
          if(is.null(input$mirexp_mouse) || 
             is.null(x <- getMouseMirExp(input$mirexp_mouse))) return(NULL)
          x <- head(x,round((input$mir_cut_off2/100)*length(x)))
          x <- setNames(rep(x,each=3),
                        paste0(rep(names(x),each=3), c("","-5p","-3p")))
          names(x) <- gsub("mir","miR",names(x))
          return(data.frame(row.names=names(x), expression=as.numeric(x)))
        }
        return(NULL)
      }
      if(is.null(input$exp_mirna_file)) return(NULL)
      mirup <- as.data.frame(data.table::fread(input$exp_mirna_file$datapath))
      # if(ncol(mirup)!=2 || !is.numeric(mirup[,2])){
      #   showModal("The miRNA data you entered is not valid")
      # }
      mirup <- mirup[order(mirup[[2]], decreasing=TRUE),]
      colnames(mirup)[1] <- "name"
      colnames(mirup)[2] <- "expression"
      mirup <- mirup[1:((input$mir_cut_off/100)*nrow(mirup)),]
      return(data.frame(row.names=mirup[,1], expression=mirup[,2]))
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
    
    observeEvent(input$example_GOI, {
      goi <- paste(.exampleGeneset(), collapse=", ")
      bg <- paste(.exampleBackground(), collapse=", ")
      if(tolower(input$species) != "human"){
        goi <- tools::toTitleCase(tolower(goi))
        bg <- tools::toTitleCase(tolower(bg))
      }
      updateTextAreaInput(session, "background_genes", value=bg)
      updateTextAreaInput(session, "genes_of_interest", value=goi)
      updateSelectInput(session, "genes_format", "GS")
    })
    
    ##############################
    ## Initialize Genes of Interest
    
    Gene_Subset <- reactive({
      if(isTRUE(getOption("shiny.testmode"))) print("Gene_Subset")
      if(is.null(input$input_type) || is.null(EN_Object())) return(NULL)
      if(input$input_type == "dea"){
        if(is.null(input$dea_input)) return(NULL)
        d <- DEA()
        return(row.names(d)[d$FDR<input$binary_fdr])
      }else{
        if(input$GOI == "GOI_custom"){
          g <- trimInputList(input$genes_of_interest)
          g <- gsub("\\..*","",g)
        }else{
          Sp <- switch(input$species, Human="Hs", Mouse="Mm", Rat="Rn")
          ens <- switch(input$genes_format, "Ens"=TRUE, "GS"=FALSE)
          return(as.character(unlist(getGOgenes(go_ids=input$go_term, 
                                                species=Sp,ensembl_ids=ens))))
        }
      }
    })
    

  
    ##############################
    ## Initialize target predictions

    
    EN_Object <- reactive({
      if(isTRUE(getOption("shiny.testmode"))) print("EN_Object")
      if(is.null(input$collection) || is.null(input$species) ||
         is.null(bDataLoaded[[input$species]][[input$collection]])) return(NULL)
      if(is.character(bDataLoaded[[input$species]][[input$collection]])){
        #waiter::waiter_show(html=tags$h3("Loading database..."), color="#333E4896")
        bDataLoaded[[input$species]][[input$collection]] <<-
          readRDS(bDataLoaded[[input$species]][[input$collection]])
        #waiter::waiter_hide()
      }
      bDataLoaded[[input$species]][[input$collection]]
    })
    
    output$collection_details <- renderText({
      if(is.null(EN_Object())) return(NULL)
      
      paste(nrow(EN_Object()), "bindings of ", length(unique(EN_Object()$set)),
            "families on ", length(unique(EN_Object()$feature)), 
            "transcripts/genes")
    })
   
    ##############################
    ## CD plot

    flags <- reactiveValues(CDplotOn=FALSE, enrichPlotOn=FALSE)
    
    observe({
      if( input$input_type!="dea" || is.null(DEA()) ){
        hideElement("cdplot_outer")
        showElement("cdplot_na")
        flags$CDplotOn <- FALSE
      }else{
        hideElement("cdplot_na")
        showElement("cdplot_outer")
        flags$CDplotOn <- TRUE
      }
    })
    observe({
      if(flags$CDplotOn && flags$enrichPlotOn){
        plotlyObs$resume()
      }else{
        plotlyObs$suspend()
      }
    })

    CDtypeOptions <- reactive({
      if(is.null(EN_Object())) return(NULL)
      CN <- colnames(EN_Object())
      if(!any(c("sites","score","best_stype","type") %in% CN))
        return(c(Automatic="auto"))
      options <- c( Automatic="auto", "Best site type"="best_stype", 
                    Score="score", "Number of sites"="sites",
                    "Site type"="type" )
      if(!("best_stype" %in% CN) && sum(grepl("[6-8]mer",CN))>1)
        CN <- unique(c(CN,"best_stype", "type"))
      options[options %in% c("auto",CN)]
    })
    
    mirfamChoices <- reactive({
      if(is.null(EN_Object())) return(c())
      names(lvl) <- lvl <- levels(as.factor(EN_Object()$set))
      if(!is.null(m <- metadata(EN_Object())$families) &&
         !any(grepl("-miR-",EN_Object()$set, ignore.case=TRUE))){
        if(length(w <- which(lvl %in% as.character(m)))>0)
          w <- names(lvl)[w]
        if(isTRUE(getOption("shiny.testmode"))) print("mir_fam1")
        x <- sapply(split(names(m), m)[w], FUN=function(x){
          gsub("/(mir|let)", "/", gsub("/(hsa-|rno-|mmu-)", "/", 
                                       paste(x, collapse="/")), ignore.case=TRUE)
        })
        x <- setNames(names(x),as.character(x))
        lvl <- c(setdiff(lvl, x),x)
        lvl <- lvl[order(names(lvl))]
      }
      lvl
    })
    
    observe({
      if(!is.null(EN_Object())){
        updateSelectInput(session, "CD_type", choices=CDtypeOptions())
        updateSelectizeInput(session, "mir_fam", choices=mirfamChoices(),
                             server=TRUE)
      }
    })
    
    observe({
      if(!is.null(CDplot_obj())){
        if(input$CD_type %in% c("auto","best_stype","type")){
          hideElement("CD_k")
        }else{
          showElement("CD_k")
        }
      }
      })


    CDplot_obj <- reactive({
      if( input$input_type=="dea" && is.null(DEA()) ) return(NULL)
      if(is.null(input$mir_fam) || input$mir_fam=="") return(NULL)
      if(isTRUE(getOption("shiny.testmode"))) print("CDplot_obj")
      set_Name <- input$mir_fam
      if(sum(EN_Object()$set==set_Name)<5) return(FALSE)
      dea <- DEA()
      TS <- EN_Object()
      dea <- .applySynonyms(dea, TS)
      if(sum(levels(TS$feature) %in% row.names(dea))<10) return(FALSE)
      legname <- names(CDtypeOptions())[CDtypeOptions()==input$CD_type]
      if(input$CD_type=="auto")
        legname <- ifelse(length(CDtypeOptions())>1, names(CDtypeOptions())[2],
                          "Target?")
      if(isTRUE(getOption("shiny.testmode"))) print(set_Name)
      p <- tryCatch(
        CDplotWrapper(dea, TS, setName=set_Name, 
                      by=input$CD_type, k=input$CD_k) + labs(colour=legname),
        error=function(e){ e })
      if(is(p,"error")) return(FALSE)
      if(input$CDplot_xaxis==0){
        q <- quantile(dea$logFC, c(0.01, 0.99))
        p <- p + xlim(q[1],q[2])
      }else{
        p <- p + xlim(-input$CDplot_xaxis, input$CDplot_xaxis)
      }
      setTheme(p, input$CDplot_theme)
    })
              
    output$cd_plot <- renderPlot({
      p <- CDplot_obj()
      validate( need(!isFALSE(p),
             "This miRNA has an insufficient number of targets in the collection.
             (This could be because the species of the collection does not 
             match that of the input).") )
      if(!is.null(logCallsFile))
        write(paste(Sys.Date(),session$token,"CDplot"), logCallsFile,
                   append=TRUE)
      p
    })
    
    output$cd_plot_dl <- downloadHandler(
      filename={
        if(is.null(CDplot_obj()) || isFALSE(CDplot_obj())) return(NULL)
        paste0("CDplot_", input$mir_fam, ".pdf")
      },
      content = function(con){
        if(is.null(CDplot_obj()) || isFALSE(CDplot_obj())) return(NULL)
        pf <- paste0(tempfile(),".pdf")
        ggsave(pf, CDplot_obj(), device="pdf", width=7, heigh=5)
        file.copy(pf, con)
      }
    )
    
    ##############################
    ## Enrichment analysis
    
    testsAvailable <- reactive({
      if(is.null(ER())) return(c())
      nn <- setdiff(names(ER()), "regmir.bb3")
      choices <- list()
      if(any(grepl("\\.up$|\\.down$",nn)))
        choices <- list("Downregulated genes"=grep("\\.down$",nn,value=TRUE),
                        "Upregulated genes"=grep("\\.up$",nn,value=TRUE))
      if(length(tt <- setdiff(grep("overlap|regmir\\.bb", nn, value=TRUE), 
                              unlist(choices)))>0)
        choices$Binary <- tt
      if(length(tt <- setdiff(nn, unlist(choices)))>0)
        choices$Continuous <- tt
      choices <- lapply(choices, FUN=function(x) setNames(x,x))
      if(length(nn)>1) choices[["All tests"]] <- c("merged"="")
      choices
    })
    
    # Get user-friendly choices just for siteoverlap
    observe({
      updateSelectInput(session, "view_test", choices=testsAvailable())
    })
    
    ER <- eventReactive(input$enrich, {
      hideElement("resultsbox")
      if(is.null(input$input_type) || is.null(EN_Object())) return(NULL)
      if(isTRUE(getOption("shiny.testmode"))) print("ER")
      mirexp <- miRNA_exp()
      if(!is.null(mirexp) && is.null(dim(mirexp))){
        mirexp <- data.frame(row.names=mirexp, miRNA=mirexp)
      }
      if(input$input_type == "dea"){
        if(is.null(DEA())) return(NULL)
        sig <- DEA()
        bg <- NULL
        standard_tests <- c("siteoverlap","areamir")
      }else{
        if(is.null(Gene_Subset()) || is.null(Back())) return(NULL)
        sig <- Gene_Subset()
        bg <- Back()
        standard_tests <- c("siteoverlap")
      }
      tests <- c(standard_tests, input$tests2run)
      if(!any(c("sites","score") %in% colnames(EN_Object())))
        tests <- "overlap"
      tests <- intersect(tests, availableTests(EN_Object()))
      msg <- paste0("Performing enrichment analyses with the following tests: ",
                    paste(tests,collapse=", "))
      message(msg)
      if(length(sig)==0) return(NULL)
      detail <- NULL
      if(input$collection == "Targetscan all miRNA BS")
        detail <- "This will take a while..."
      if(!is.null(logCallsFile))
        write(paste(Sys.Date(),session$token,"enrich"), logCallsFile,
                   append=TRUE)
      res <- withProgress(message=msg, detail=detail, value=1, max=3, {
        tryCatch(testEnrichment(sig, EN_Object(), background=bg, 
                                sets.properties=mirexp, tests=tests, 
                                minSize=input$minsize, th.FDR=input$dea_sig_th),
                error=function(e){
                  showModal(modalDialog(
                    title = "There was an error with your request:",
                    tags$pre(as.character(e)),
                    "This typically happens when there is a mismatch 
                    between your different input data (e.g. wrong species or the first column of a DEA does not contain genes)"
                  ))
                })
      })
      showElement("resultsbox")
      res
    })
    
    #adapt plot_size
    jqui_resizable(ui="#bubble_plot")
    
    erRes <- reactive({
      if(is.null(ER())) return(NULL)
      if(is.null(test <- input$view_test) || test==""){
        er <- getResults(ER(), getFeatures=FALSE, flatten=TRUE)
        er$FDR <- er$FDR.geomean
        er$enrichment <- rowMeans(er[,grep("nrichment|beta|coefficient",
                                           colnames(er)),drop=FALSE],na.rm=TRUE)
      }else{
        er <- getResults(ER(), test=test, getFeatures=FALSE, flatten=TRUE)
      }
      er      
    })
    
    output$bubble_plot <- renderPlotly({
      if(is.null(er <- erRes())) return(NULL)
      if(isTRUE(getOption("shiny.testmode"))) print("bubble_plot")
      col.field <- "expression"
      if(is.null(er$expression)) col.field <- NULL
      er$set <- row.names(er)
      flags$enrichPlotOn <- TRUE
      p <- enrichPlot(er, repel=FALSE, label.sig.thres=input$label.sig.thres,
                      sig.field=input$sig.field, col.field=col.field, 
                      label.enr.thres=input$label.enr.thres, 
                      maxLabels=input$label_n )
      p <- setTheme(p, input$bubble_theme)
      forTooltip <- intersect(c("set","label","overlap","enrichment","set_size",
                                "pvalue","FDR"), colnames(er))
      ggplotly(p, source="enrichplot", tooltip=forTooltip)
    })
    
    plotlyObs <- observeEvent(event_data("plotly_click", "enrichplot",
                                          priority="event"), suspended=TRUE, {
      if(is.null(er <- erRes()) || !flags$CDplotOn || is.null(input$mir_fam))
        return(NULL)
      event <- event_data("plotly_click", "enrichplot")
      if(!is.list(event) || is.null(event$pointNumber)) return(NULL)
      rid <- as.integer(event$pointNumber+1)
      if(is.null(rid) || !(rid>0)) return(NULL)
      fam <- row.names(er)[rid]
      if(!(fam %in% mirfamChoices()))
        fam <- strsplit(er[rid,"members"],";")[[1]][[1]]
      if(!(fam %in% mirfamChoices())) return(NULL)
      updateTabItems(session, "main_tabs", "tab_cdplot")
      updateSelectizeInput(session, "mir_fam", choices=mirfamChoices(), 
                           server=TRUE, selected=fam)
    })
    
    output$bubble_plot_dl <- downloadHandler(
      filename={
        if(is.null(erRes())) return(NULL)
        paste0("enrichPlot_", input$view_test, ".pdf")
      },
      content = function(con){
        if(is.null(er <- erRes())) return(NULL)
        col.field <- "expression"
        if(is.null(er$expression)) col.field <- NULL
        er$set <- row.names(er)
        
        p <- enrichPlot(er, repel=TRUE, label.sig.thres=input$label.sig.thres,
                        sig.field=input$sig.field, col.field=col.field, 
                        label.enr.thres=input$label.enr.thres, 
                        maxLabels=input$label_n )
        p <- setTheme(p, input$bubble_theme)
        pf <- paste0(tempfile(),".pdf")
        ggsave(pf, p, device="pdf", width=8, heigh=5)
        file.copy(pf, con)
      }
    )
    
    output$hoverinfo <- renderUI({
      id <- suppressWarnings(tryCatch(
                        event_data(event="plotly_hover", source="enrichplot"),
                        error=function(e) return(tags$p(HTML("&nbsp;")))))
      if(is.null(id) || is.null(ER())) return(tags$p(HTML("&nbsp;")))
      test <- input$view_test
      if(is.null(test) || test=="") test <- NULL
      rr <- getResults(ER(), test=test, flatten=TRUE, getFeatures=FALSE)
      rr <- rr[id[1,2]+1,]
      if(is.null(rr$members)) return(tags$p(row.names(rr)))
      tags$p(tags$b("Members: "), gsub(";", "; ", rr$members))
    })
    
    
    output$hits_table <- renderDT({ # prints the current hits
      if(is.null(ER())) return(NULL)
      test <- input$view_test
      if(is.null(test) || test=="") test <- NULL
      rr <- getResults(ER(), test=test, flatten=TRUE)
      show_standard <- c("enrichment","pvalue","FDR","expression")
      if(is.null(input$columns2show)){
        columns2hide <- setdiff(colnames(rr),show_standard)
      }else{
        show_add <- input$columns2show
          if(any(show_add %in% "genes")){
            show_add <- c(show_add,"genes.down","genes.up","genes.features")
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
        fn <- paste0("EnrichMir_hits_",input$view_test,"_",Sys.Date(),".csv")
        fn
      },
      content = function(con) {
        if(is.null(ER())) return(NULL)
        test <- input$view_test
        if(is.null(test) || test=="") test <- NULL
        h <- getResults(ER(), test=test, flatten=TRUE)
        write.csv(h, con)
      }
    )
      
    waiter::waiter_hide()
    if(isTRUE(getOption("shiny.testmode"))) print("END LOADING")
    
  }
}

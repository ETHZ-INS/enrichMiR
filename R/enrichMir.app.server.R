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
                            columnDefs=list(
                              list(visible=FALSE, 
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
      ##if (is.null(upFile)) return(NULL)
      if (is.null(upFile)){ ## TEMPORARY - BECAUSE I'M LAZY
        updf <- readRDS("/mnt/schratt/enrichMiR_benchmark/data/bartel_HEK.deaList.rds")[[1]]
      }else{
        updf <- data.table::fread(upFile$datapath)
      }
      updf <- enrichMiR:::.homogenizeDEA(updf)
      if(nrow(updf)==0) return(NULL)
      updf
    })
    
    output$dea_res <- renderUI({
      if(is.null(DEA()))
        return(tags$span(icon("exclamation-triangle"), "No valid file uploaded",
                         style="font-weight: bold; font-color: red;"))
      tags$span(icon("check-circle"), "Valid file",
                style="font-weight: bold; font-color: forestgreen;")
    })
    
    output$enrichbtn <- renderUI({
      if( (input$input_type == "dea" && is.null(DEA())) ||
          (input$input_type != "dea" && (
            is.null(Gene_Subset()) || length(Gene_Subset())<2 ||
            is.null(Back()))) )
        return(tags$span(icon("exclamation-triangle"), "No valid gene/DEA input!",
                         style="font-weight: bold; font-color: red;"))
      actionButton(inputId = "enrich", "Enrich!", icon = icon("search"))
    })
    
    ## Include , and "" gsub
    Back <- reactive({ #initalize the background
      if(input$input_type == "dea"){
        if(is.null(DEA())) return(NULL)
        row.names(DEA())
      }else{
        b <- trimInputList(input$background_genes)
        b <- gsub("\\..*","",b)
      }
    })
        
    miRNA_exp <- reactive({
      if(isTRUE(getOption("shiny.testmode"))) print("miRNA_exp")
      if(is.null(input$expressed_mirna_box)) return(NULL)
      if(input$expressed_mirna_box=="Custom Set"){
        if(is.null(input$exp_mirna_list) || input$exp_mirna_list=="") return(NULL)
        return(trimInputList(input$exp_mirna_list))
      }
      if(is.null(input$exp_mirna_file)) return(NULL)
      mirup <- as.data.frame(data.table::fread(input$exp_mirna_file$datapath))
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
    
    observeEvent(input$example_GOI, {
      goi <- "PLEKHB2,FBXO21,PIGS,SLC7A1,YKT6,RABL6,SLC1A5,OCRL,SH3RF1,SLC9A1,SLC7A5,SRPRA,G6PC3,ARL2,CTDSPL,CSRP1,TBC1D22B,ZER1,SLC52A2,GNPDA1,PRDM4,VAMP3,TMEM87A,SEPTIN2,SHCBP1,RELL1,PRTFDC1,CPEB1,PHKA1,RNF38,TMED8,CASP7,PROSER1,LMNB2,DSTYK,TMEM216,RAVER1,PRUNE1,ASF1B,NCDN,TGFBRAP1,CS,MTHFD2,SFT2D1,ARHGAP1,IQGAP1,ATN1,CTDNEP1,VPS37B,PLCD3,PKM,POLR3D,SLC25A6,PRKRA,HECTD3,SULF1,SERINC5,DYNC1LI2,BCAT2,VPS4A"
      bg <- "PLEKHB2,FBXO21,PIGS,SLC7A1,YKT6,RABL6,SLC1A5,OCRL,SH3RF1,SLC9A1,SLC7A5,SRPRA,G6PC3,ARL2,CTDSPL,CSRP1,TBC1D22B,ZER1,SLC52A2,GNPDA1,PRDM4,VAMP3,TMEM87A,SEPTIN2,SHCBP1,RELL1,PRTFDC1,CPEB1,PHKA1,RNF38,TMED8,CASP7,PROSER1,LMNB2,DSTYK,TMEM216,RAVER1,PRUNE1,ASF1B,NCDN,TGFBRAP1,CS,MTHFD2,SFT2D1,ARHGAP1,IQGAP1,ATN1,CTDNEP1,VPS37B,PLCD3,PKM,POLR3D,SLC25A6,PRKRA,HECTD3,SULF1,SERINC5,DYNC1LI2,BCAT2,VPS4A,AL356776.1,DMAC1,MROH8,MSH2,PANX2,PAQR9,SWI5,DHX9P1,RAP1GDS1,OLFM2,PGF,TBC1D12,TNFSF12,ANP32BP1,OSGIN2,LMCD1,CNOT9,TTC6,YTHDF2,BOD1,ZMYM3,USF1,STK4,PYCR3,AC104083.1,NSUN5P2,PTAFR,CFAP58,CDKN2AIP,ARPC2,HMGN2P46,PJA1,HSPA6,AC079416.3,MYH9,MRTFB,SDK1,HNRNPA1P69,RAB10,BUD13,NDUFA2,STK31,SEC11B,ROMO1,KDELC1P1,EXOSC5,AC092807.3,AC004492.1,PRSS53,MTRNR2L6,RPS2P6,RPL5P1,AC013391.2,PNPT1P1,CALM2,AC027801.2,ANO8,RAB2B,AL136038.4,PPM1B,BHLHE41,SBF1,LIFR-AS1,PSMB4,POGK,CSDC2,PPIAP26,PHRF1,NPM1P9,PDIA2,BCRP2,GAPDHP58,CLEC16A,AC234775.3,EP400,VTI1B,PTPN6,SIMC1,AL034397.1,RBBP8,PEX14,ABCA11P,TEN1-CDK3,ZYG11B,AGO1,TAF3,GBX2,PLIN5,PKDCC,SMNDC1,RUSC2,FAM20C,HBP1,SDHAF1,ZBTB40-IT1,TMEM14C,TUFM,SH3GL2,MDM4,MIR3936HG,FAM210A,SNRPGP10,CCAR1,RORB,AL645940.1,RPN1,LSM6,PPRC1,CD70,AC027307.3,ZNF789,BRSK2,MEGF6,DPYSL3,CRYBG3,DMAP1,TMEM101,CDKN1B,UCP1,ZNF222,LINC00926,RGS7,EML4,PCBP1,BTF3P8,INTS1,AC106786.1,RNF168,CASP9,MLXIPL,CDX2,ANKRD54,AC027682.7,PRRG1,AC209007.1,AC026356.1,AC104964.3,THG1L,AC139530.1,TEX22,RPL7P59,CRB2,PDCD2,HTR1D,PPIA,YARS2,AC092183.1,TVP23C,EEF1A1P6,DNAAF1,INAFM1,TRIOBP,PIR,NOP10,DDX50P2,ACTL8,UQCRC2,RAB7A,FBXL19-AS1,PA2G4P4,MAP4K2,HADHAP2,CKAP5,ZNF215,GPX4,MTX1,TTLL4,AHNAK,KRT18P31,ELP4,NCAM1,DOK4,ACAA2,AC092431.2,RPL22P8,SP9,ST6GALNAC6,TRA2A,NOP2,PLXNB2,GADD45GIP1,ASH1L-AS1,ZNF843,RPS23P8,KARS,ZNF407,TRAPPC2B,ZNF518A,KCNA3,UCKL1,NPM1P6,NOL3,SPR,NOP14,LGMN,VWA5A,GAPDHP25,ARTN,KBTBD7,DCXR,ANO6,AC008147.4,HAX1,NEDD1,ZNF772,AC020910.5,DEPDC7,HERC2P2,DCHS1,RPL37P1,PLIN1,APLF,ASB13,NOP56,PCSK4,AMACR,PLCL1,G0S2,VASN,PPP3CB,PSMD4,AL512791.2,STK17A,TSSK3,AC006460.2,PPP1R26,ELF1,ACTBP8,CSNK2B,EDRF1,AC027097.1,PTBP1,KLLN,GSTK1,VKORC1L1,COL11A2,C15orf62,PRKCG,ACVR2B-AS1,RDH16,WASHC2C,AC114980.1,TXN2,RPL19P16,SRR,RPS19BP1,AC010680.5,RN7SL832P,PCNA,AC013403.2,UBE2Q2,AC016737.1,TMPO,API5P2,FAAH2,STAMBP,AL445487.1,RPS2P44,SRRM2-AS1,AC079193.2,GATAD2A,FEZF1,HHEX,RPF1,IRX3,AC093788.1,STON2,TCN2,CMSS1,RPL27AP5,DPH1,AL136116.3,TOMM5,RPL14,AL365295.1,AC090015.1,VPS13B,MXD4,CHP2,RN7SL535P,ACTG1P14,AC087343.1,PRPF38AP2,GRIPAP1,AC107068.1,CTSL,CTSF,SNHG11,HHIP-AS1,NMUR1,R3HCC1L,PPP1R13L,UGDH-AS1,PDF,GINS4,PIPSL,AC109347.1,SLC16A10,RPS18P13,ASPHD2,EEF1E1P1,AL096803.2,RIC8B,TFDP1,RBM6,METTL4,SULT4A1,RPS2P46,DSTN,AC005072.1,ACAA1,AC012313.1,GSDME,HNRNPA1P54,TOGARAM2,AC022498.2,TVP23B,RPS12P26,AC131235.1,RCCD1,INPP1,WDFY1,AC026367.2,MCF2L-AS1,SETP21,TMEM138,SNORA33,PABPN1,OGDH,AL122023.1,CHEK1,CCN2,AL139260.1,SYMPK,CDKN2AIPNL,MBD5,NCOA4P4,BLVRB,DSE,MIEF2,RBM8B,ZNF350,ADSL,PARN,LNCTAM34A,PNMA6A,PIGW,TIMM29,IL21R-AS1,NBPF3,CYB5RL,CERT1,MAP7D2,RUNX1,AC127024.8,JAKMIP2,NIBAN2,TIGAR,RNF227,UNG,MTOR,TMEM107,JUND,GEMIN7,CALCOCO2,TOR2A,ATAD3A,C1orf43,CRK,AL353151.1,GAPDHP23,DEGS1,KXD1,AF131216.4,GTF2H2C,LAMP2,IPO7P1,RAP2C-AS1,LFNG,MTND2P9,CNOT7,AC134349.1,GMPSP1,NAA25,MMP24OS,VMP1,NTNG2,MCCC1,RPS25,ILF2P2,AL035530.2,AC024451.1,GNMT,XRCC6P1,NSMCE3,AC026124.1,AC074194.1,NHP2P1,CHD2,AC073046.1"
      updateTextAreaInput(session, "background_genes", value=bg)
      updateTextAreaInput(session, "genes_of_interest", value=goi)
      updateSelectInput(session, "genes_format", "GS")
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
          g <- trimInputList(input$genes_of_interest)
          g <- gsub("\\..*","",g)
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
                                                   "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds"),
                                                   "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_mouse.rds"),
                                                   "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_rat.rds"),
                                                   "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds")),
             "Targetscan conserved miRNA BS" = switch(input$species,
                                                        "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds"),
                                                        "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_mouse.rds"),
                                                        "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_rat.rds"),
                                                        "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds")),
             "Targetscan all miRNA BS" = switch(input$species,
                                                       "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_AllSites_human.rds"),
                                                       "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_AllSites_mouse.rds"),
                                                       "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_AllSites_rat.rds"),
                                                       "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_AllSites_human.rds")),
             "CISBP RBP motif sites" = switch(input$species,
                                                         "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds"),
                                                         "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_mouse.rds"),
                                                         "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_rat.rds"),
                                                         "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds")),
             "miRTarBase" = switch(input$species,
                                              "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds"),
                                              "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_mouse.rds"),
                                              "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_rat.rds"),
                                              "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds")),
             "Custom - not yet" =  switch(input$species,
                                          "Human" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds"),
                                          "Mouse" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_mouse.rds"),
                                          "Rat" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Mouse_ConSites_rat.rds"),
                                          "Custom - not yet" = readRDS("/mnt/schratt/enrichMiR_data/Targetscan/20201102_Targetscan_Human_ConSites_human.rds")))
    })

                                 
   
    
   
    ##############################
    ## CD plot
    
    observe({
      if(!is.null(EN_Object())){
        if(!is.null(m <- metadata(EN_Object())$families)){
          if(!is.null(miRNA_exp())) {
            mirexp <- miRNA_exp()
              if(is.null(dim(mirexp))) {
                mirexp <- mirexp
              }else{
                mirexp <- mirexp[,1]
              }
            m <- m[names(m) %in% mirexp]
          }
        updateSelectizeInput(session, "mir_fam", choices=m, server=TRUE)
        }
      }
    }) 
              
    output$cd_plot <- renderPlot({
      if(is.null(input$mir_fam) || input$mir_fam=="") return(NULL)
      validate(
        need(any(EN_Object()$set %in% input$mir_fam), "This miRNA has no annotated Binding Site")
      )
      dea <- DEA()
      TS <- EN_Object()
      dea <- .applySynonyms(dea, TS)
      CDplot2(dea, TS, setName=input$mir_fam, by = input$CD_type, k=input$CD_k) + 
        xlim(-input$CDplot_xaxis, input$CDplot_xaxis)
    })
    
    ##############################
    ## Enrichment analysis
    
    # Get choices for all tests
    observe({
      if(!is.null(ER())) updateSelectInput(session, "view_test", choices=c("",names(ER())))
    })
    
    testsAvailable <- reactive({
      if(is.null(ER())) return(c())
      nn <- names(ER())
      choices <- list()
      if(any(grepl("\\.up$|\\.down$",nn)))
        choices <- list("Downregulated genes"=grep("\\.down$",nn,value=TRUE),
                        "Upregulated genes"=grep("\\.up$",nn,value=TRUE))
      if(length(tt <- setdiff(grep("overlap|regmir", nn, value=TRUE), unlist(choices)))>0)
        choices$Binary <- tt
      if(length(tt <- setdiff(nn, unlist(choices)))>0)
        choices$Continuous <- tt
      choices <- lapply(choices, FUN=function(x){ names(x) <- x; x})
      if(length(nn)>1) choices[["All tests"]] <- c("merged"="")
      if(isTRUE(input$view_all)) return(choices)
      choices <- lapply(choices, FUN=function(x) 
        head(grep("siteoverlap|areamir",as.character(x),value=TRUE),1))
      choices[lengths(choices)>0]
    })
    
    # Get user-friendly choices just for siteoverlap
    observe({
      updateSelectInput(session, "view_test", choices=testsAvailable())
    })
    
# 
#     output$test_info <- renderPrint({ # print the test
#       if(is.null(ViewTest()) || ViewTest()=="") return("")
#       out <- capture.output(ViewTest())
#       cat(out)
#     })
#       
    
    ER <- eventReactive(input$enrich, {
      if(is.null(input$input_type) || is.null(EN_Object())) return(NULL)
      if(isTRUE(getOption("shiny.testmode"))) print("ER")
      mirexp <- miRNA_exp()
      if(!is.null(mirexp) && is.null(dim(mirexp)))
        mirexp <- data.frame(row.names=mirexp, members=mirexp)
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
      msg <- paste0("Performing enrichment analyses with the following tests: ",paste(tests,collapse=", "))
      message(msg)
      if(length(sig)==0) return(NULL)
      detail <- NULL
      if(input$collection == "Targetscan all miRNA BS") detail <- "This will take a while..."
      withProgress(message=msg, detail=detail, value=1, max=3, {
        o <- tryCatch(testEnrichment(sig, EN_Object(), background=bg, 
                                     sets.properties=mirexp, tests=tests, 
                                     minSize=input$minsize, th.FDR=input$dea_sig_th),
                      error=function(e){
                        showModal(modalDialog(
                          title = "There was an error with your request:",
                          tags$pre(as.character(e)),
                          "This typically happens when there is a mismatch between your different input data (e.g. wrong species)"
                        ))
                      })
      })
    })
    
    output$bubble_plot <- renderPlotly({
      if(is.null(ER())) return(NULL)
      test <- input$view_test
      if(is.null(test) || test==""){
        er <- getResults(ER(), getFeatures=FALSE, flatten=TRUE)
        er$FDR <- er$FDR.geomean
        er$enrichment <- rowMeans(er[,grep("nrichment|beta|coefficient",colnames(er)),drop=FALSE],na.rm=TRUE)
      }else{
        er <- getResults(ER(), test=test, getFeatures=FALSE, flatten=TRUE)
      }
      ggplotly(enrichPlot(er, repel=FALSE, label.sig.thres = input$label.sig.thres, sig.field = input$sig.field,
                          label.enr.thres = input$label.enr.thres, maxLabels = input$label_n ))
    })
    
    
    output$hits_table <- renderDT({ # prints the current hits
      if(is.null(ER())) return(NULL)
      test <- input$view_test
      if(is.null(test) || test=="") test <- NULL
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
      
    
    if(isTRUE(getOption("shiny.testmode"))) print("END LOADING")
    
  }
}





  
  
  



#' @import shiny DT shinydashboard shinycssloaders
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyjqui jqui_resizable
#' @importFrom shinyjs useShinyjs
#' @importFrom rintrojs introjsUI
#' @importFrom waiter use_waiter waiter_show_on_load waiter_hide
#' @export
enrichMiR.ui <- function(){
  library(shiny)
  library(shinyjqui)
  library(DT)
  library(shinydashboard)
  library(shinycssloaders)
  library(plotly)
  library(ggplot2)
  library(rintrojs)
  library(waiter)
  library(shinyjs)

  genes_placeholder <- paste("Enter your genes as symbols or ensembl IDs,",
                             "separated by spaces, commas, or linebreaks. E.g.:
EZH2, YY1, SHANK3, ...\nor:
ENSG00000106462, ENSG00000100811, ...")
  ggplot_themes <- setdiff(grep("^theme_",ls(getNamespace("ggplot2"), 
                                             all.names=TRUE),value=TRUE), 
                           paste0("theme_",c("all_null","set","update",
                                             "get","replace","void")))
    
  ui <- dashboardPage(
    dashboardHeader(title = "enrichMiR", titleWidth = "300px",
                    tags$li(class="dropdown",
                            actionLink("helpBtn", label="Help",
                                       icon=icon("question")))),
    ## Sidebar content
    dashboardSidebar(width = "300px",
       sidebarMenu(
         menuItem("Introduction", tabName = "tab_intro", icon=icon("info")),
         tags$hr(width="80%"),
         menuItemOutput("menu_species"),
         menuItemOutput("menu_input"),
         menuItemOutput("menu_enrich"),
         menuItemOutput("menu_cdplot"),
         tags$hr(width="80%"),
         menuItem("Tests benchmark", tabName="tab_benchmark",
                  icon=icon("tachometer-alt"))
       )                     
    ),
    
    ## Body Content
    dashboardBody(
#       tags$head(tags$style(HTML("
# #sel_test_div label {
#   display: table-cell;
#   text-align: center;
#   vertical-align: middle;
# }
# #sel_test_div .form-group {
#   display: table-row;
# }
# "))),
      useShinyjs(), introjsUI(), use_waiter(spinners = 3),
      waiter_show_on_load(html=tagList(
        #tags$h1("enrichMiR"),
        tags$img(src="http://130.60.24.189:81/common/enrichMiR/enrichMiR_sticker.png"),
        tags$h3("Please wait while the application is initialized..."),
        spin_1()
      )),
      tabItems(
        tabItem(tabName = "tab_intro",
          box(title="enrichMiR: miRNA (and RBP) target enrichment analysis", 
              width=12, tags$p(
                "This app will allow you to identify miRNAs whose targets are ",
                "enriched among genesets of interest or a differential ",
                "expression signature, and produce related visualizations.",
                "Although the app was chiefly developed (and benchmarked) for ",
                "miRNAs, some support is also offered to run the same analyses",
                " using RNA-binding proteins."), tags$br(),
              tags$p("To get started, take a ", 
                     actionLink("helpLink", "quick tour"),
                     " of the app."), tags$br(),
              tags$p( style="text-align: right;",
                      paste("enrichMiR version",packageVersion("enrichMiR")),
                      "; ", tags$a( href="http://schrattlab.ethz.ch",
                                    "Schratt lab", target="_blank") )
          )
        ),
        
        tabItem(tabName = "tab_species", ############ SPECIES / COLLECTION
          box(title="Select Species and Collection", width=12,
            column(7, selectInput(inputId = "species", "Species", width='98%',
                        choices = c("Human", "Mouse", "Rat")), tags$br(),
               tags$div(id="collection_input",
                  actionButton(inputId="help_collections", style="float:right;",
                           icon=icon("question-circle"), label=""),
                  selectInput(inputId="collection", width='98%', choices=c(),
                    label="Select a binding sites collection",
                    selected="Targetscan conserved miRNA BS")
                  )),
            column(5, withSpinner(textOutput("collection_details")))),
          tags$div(id="exprMirs_box", box(title="Specify expressed miRNAs", 
              collapsible=TRUE, collapsed=TRUE, width=12,
            tags$p("miRNA expression can be used to restrict and annotate the ",
                   "enrichment analysis. This information can either be given ",
                   "manually, provided by uploading a miRNA profile, or ",
                   "fetched from pre-loaded miRNA expression profiles"),
            tags$div(id="expressed_mirna_outer", 
            tabBox(id="expressed_mirna_box", width=12,
              tabPanel(title = "Custom Set",
                       "Paste a list of expressed miRNAs in 'miRBase' format", br(), br(),
                        textAreaInput(inputId = "exp_mirna_list", 
                            label="miRNA List", 
                            rows=5,
                            placeholder="Enter miRbase IDs, separated by spaces, commas, or linebreaks. E.g.:
  hsa-miR-30b-5p
  hsa-miR-30d-5p
  ...", 
                            resize="vertical")),
              tabPanel(title = "Upload miRNA expression table",
                "Upload miRNA Expression Table in the following format: miRBase name in the first
                            column, expression-values in the second column", br(),
                fileInput(inputId="exp_mirna_file", 
                          label="Upload miRNA expression object as '*.csv' file (see help)",
                          accept=c("text/csv", ".csv",".tab",".txt",
                                   "text/comma-separated-values,text/plain")),
                sliderInput(inputId="mir_cut_off", 
                            label="Select miRNA expression cut-off:",
                            min=10, max=100, post=" %", value=50),
                "Keep only the top ..% expressed miRNAs"
              ),
              tabPanel(title="Use preset expression profile",
                withSpinner(uiOutput("mirexp_preset")),
                sliderInput(inputId="mir_cut_off2", 
                            label="Select miRNA expression cut-off:",
                            min=10, max=100, post=" %", value=50 ),
                "Keep only the top ..% expressed miRNAs"
              )
            )), br(),
            footer="Note: If no miRNAs are uploaded/selected, enrichment 
                    searches will be performed with all microRNAs of the given species"))
        # eventually enable this as upload or selection from pre-loaded tissues?
        ),
        
        
        tabItem("tab_input", ############ INPUT
          tabBox(id="input_type", width=12,
            tabPanel(title = "Select geneset & background", 
              tags$p("In this mode, your genes of interest are compared against a background of genes."),
              tags$h3("Genes of interest"),
              tabsetPanel(id="GOI",
                tabPanel(title = "Custom set", value = "GOI_custom",
                         "Paste a list of expressed genes in the selected format",
                         br(), br(),
                         textAreaInput(inputId = "genes_of_interest", 
                                       label="Gene List", rows=5, width="100%",
                                       value=NULL, resize="vertical",
                                       placeholder=genes_placeholder),
                ),
                tabPanel(title = "From Gene Ontology", value="GOI_GO",
                         radioButtons(inputId = "genes_format", label = "Select:", 
                                      choices = c("Ensembl" = "Ens",
                                                  "Gene Symbol" = "GS"),
                                      selected = "Ens"),
                         selectizeInput(inputId="go_term", "Select GO-Term", choices=c()),
                         textOutput("GOI_nb"))
              ),
              br(), tags$hr(), tags$h3("Background"),
              textAreaInput(inputId = "background_genes", 
                            label="Gene List", 
                            rows=5, width="100%",
                            value=NULL,
                            placeholder=paste(genes_placeholder,"
Note that the background should also include the genes of interest!"),
                            resize="vertical"),br(),
              actionButton("example_GOI", "Example genes"),
              footer=paste("Note: If you want to use the Targetscan miRNA annotations",
                           " together with rat genes, use the 'Gene Symbol' format")
            ),
            tabPanel(title = "Upload DEA results", value = "dea",
              fluidRow(       
                column(8, fileInput( inputId = "dea_input", 
                  label="Upload DEA object as '*.csv' file (see help)",
                  accept=c("text/csv", ".csv", ".tab", ".txt",
                           "text/comma-separated-values,text/plain")),
                       "Or: ", actionButton("example_dea", "Use example DEA"),
                  tags$hr(), tags$br(),
                  tags$p("Upload Differential Expression Analyses (DEAs) as table
                with at least following information: Provide ENSEMBL_ID or Gene 
                Symbol as identifier in the first column (same format as for 
                the background), as well as logFC-values and FDR-values")),
                column(4, htmlOutput("dea_res"), tags$br(), tags$hr(), tags$br(),
                     sliderInput(inputId = "dea_sig_th", 
                                 label="Select significance threshold",
                                 min=0.01, max=0.5, value=0.05, step=0.01))
              )
            )
          )
        ),
        
        
        tabItem(tabName = "tab_enrich", ###### ENRICHMENT ANALYSIS
          column(3, tags$div(id="enrichbtn_outer", 
                             withSpinner(uiOutput("enrichbtn")))), 
          box(width=9, title="Enrichment options", collapsible=TRUE, collapsed=TRUE,
              sliderInput(inputId="minsize", 
                          label=paste(
                            "Minium number of annotated targets that is required to",
                            "consider the miRNA-family for testing"),
                          min=1, max=20, value=5, step=1),
              br(),
              box(width="100%", title = "Advanced Options", collapsible=TRUE, collapsed=TRUE, 
                  tags$h5(em(
                    "We recommend to only change test settings after reading the",
                    " enrichMir documentation and the benchmark.")),
                  tags$h5(em("By default, we always perform the 'siteoverlap'",
                             "and (with DEA inputs) the 'areamir' tests")),
                  checkboxGroupInput(inputId="tests2run",
                                     label="Select additional tests to run",
                                     choices=list(
                                       "overlap" = "overlap",
                                       "weight.overlap" ="woverlap",
                                       "plMod" = "plMod",
                                       "modscore" = "modscore",
                                       "ks" = "ks",
                                       "mw" = "mw",
                                       "regmir" = "regmir"),
                                     selected=NULL, inline=FALSE, width=NULL)
              )
            ),
          tags$div(id="resultsbox", style="display: none; height: 100%;", 
          box(width=12, title="Results", collapsible=TRUE,
            fluidRow(
              column(7, id="sel_test_div", 
                     selectInput("view_test", "Select test to visualize", 
                                 choices=c())),
              column(1, style="margin-top: 25px;", 
                     actionButton(inputId="help_tests", 
                                  icon=icon("question-circle"), label="")),
            ),
            tabBox(id="enrichresults_out", width=12,
              tabPanel(title="Enrichment plot", 
                fluidRow(
                  column(11, tags$p("Hover on a point to view family members ",
                                    " and enrichment-related statistics.")),
                  column(1, actionButton(inputId="help_enrichplot", 
                                      icon=icon("question-circle"), label=""))),
                withSpinner(jqui_resizable(plotlyOutput("bubble_plot"))), 
                tags$br(), htmlOutput("hoverinfo"), tags$br(), 
                box(title="Plot options", width=12, collapsible=TRUE, collapsed=TRUE,
            column(6,sliderInput(inputId="label.sig.thres",
                                 "Significance threshold to display labels",
                                 min=0, max=0.25, value=0.05, step=0.01)),
            column(6,sliderInput(inputId="label.enr.thres",
                                 "Enrichment threshold to display labels",
                                 min=0.5, max=10, value=1, step=0.5)),
            column(6, numericInput(inputId="label_n", value=10, min=1, max=50,
                                   label="Max number of Labels")),
            column(6, radioButtons(inputId="sig.field", 
                                   label="Display on y-axis:",
                                   choices=c("pvalue","FDR"), selected="FDR")),
            column(6, selectInput("bubble_theme", "Theme", 
                                  choices=ggplot_themes))
                )),
              tabPanel(title="Results table", 
                fluidRow(
                 column(8, checkboxGroupInput(
                   inputId="columns2show", selected=NULL, inline=TRUE,
                   label="Select add. columns to be shown:",
                   choices=list("miRNA names"="members",
                                "predicted target genes"="genes") )),
                 column(4, downloadLink('dl_hits', label = "Download all"))),
                 withSpinner(DTOutput("hits_table")), tags$br()
              )
            )
          ))
        ),
        
        
        tabItem(tabName = "tab_cdplot",  ######### CD plot
          box(width=12, title="CD Plot", 
              tags$h3(id="cdplot_na", style="font-color: red;", 
                      icon("exclamation-triangle"), 
                      "Cumulative distribution plots require the use of a DEA input."),
              tags$div(id="cdplot_outer", style="display: none;",
                fluidRow(
                  column(8,selectizeInput(inputId="mir_fam", choices=c(),
                                      label="Select miRNA family to display")),
                  column(3,selectInput(inputId="CD_type", "Split by",
                                   choices=c("sites","score"))),
                  column(1, actionButton(inputId="help_cdplot", 
                                    icon=icon("question-circle"), label=""))),
                withSpinner(jqui_resizable(plotOutput("cd_plot", width='100%', 
                                                      height='400px'))),
                tags$br(), tags$br(), tags$br(), tags$br(), 
                box(width=12, title="Plotting options", collapsible=TRUE, 
                    collapsed=TRUE,
                  column(6,sliderInput(inputId="CDplot_xaxis",
                                       "logFC to display on x.axis",
                                       min=0, max=8, value=0, step=0.5)),
                  column(6,sliderInput(inputId="CD_k", "Approximate number of sets",
                                       min=2, max=6, value=3, step=1)),
                  column(6, selectInput("CDplot_theme", "Theme", 
                                        choices=ggplot_themes)),
                  column(6, actionButton("cd_plot_dl", "Download"))                  
                )
              )
          )
        ),
        tabItem(tabName = "tab_benchmark", "Forthcoming")
      )
    )
  )
}           




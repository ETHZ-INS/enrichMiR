#' @import shiny shinydashboard shinycssloaders
#' @importFrom plotly plotlyOutput
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyjs useShinyjs runjs
#' @importFrom shinyjqui jqui_resizable
#' @importFrom rintrojs introjsUI
#' @importFrom waiter use_waiter waiter_show_on_load waiter_hide spin_1
#' @importFrom DT DTOutput
#' @export
enrichMiR.ui <- function(){

  genes_placeholder <- paste("Enter your genes as symbols or ensembl IDs,",
                             "separated by spaces, commas, linebreaks or directly paste a spreadsheet column. E.g.:
EZH2, YY1, SHANK3, ...\nor:
ENSG00000106462, ENSG00000100811, ...")
  ggplot_themes <- c(gsub("^theme_","",setdiff(grep("^theme_",ls(getNamespace("ggplot2"), 
                                             all.names=TRUE),value=TRUE), 
                           paste0("theme_",c("all_null","set","update","gray",
                                             "get","replace","void")))),
                     "classic2", "pubr", "pubclean")
    
  ui <- dashboardPage(
    dashboardHeader(title="enrichMiR", titleWidth="300px",
      dropdownMenu(icon=tagList(icon("question"),tags$span(" Help")), 
        headerText="Documentation topics:", type="notifications", badgeStatus=NULL,
        tags$li(tags$a("Interactive tour", id="helpBtn", href="#", 
                       class="action-button shiny-bound-input")),
        tags$li(tags$a("Overview", id="overview", href="#", 
                       class="action-button shiny-bound-input")),
        tags$li(tags$a("Binding sites collections", id="help_collections2", href="#", 
                       class="action-button shiny-bound-input")),
        tags$li(tags$a("Enrichment tests", id="help_tests2", href="#", 
                       class="action-button shiny-bound-input")),
        tags$li(tags$a("CD plots", id="help_cdplot2", href="#", 
                       class="action-button shiny-bound-input")),
        tags$li(tags$a("Browser compatibility", id="brow_comp", href="#", 
                       class="action-button shiny-bound-input"))
        )
    ),
    ## Sidebar content
    dashboardSidebar(width = "300px",
       sidebarMenu(id="main_tabs",
         menuItem("Introduction", tabName = "tab_intro", icon=icon("info")),
         tags$hr(width="80%"),
         menuItemOutput("menu_species"),
         menuItemOutput("menu_input"),
         menuItemOutput("menu_enrich"),
         menuItemOutput("menu_cdplot"),
         tags$hr(width="80%"),
         menuItem("Tests & benchmark", tabName="tab_benchmark",
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
      useShinyjs(), introjsUI(), use_waiter(),
      
      waiter_show_on_load(html=tagList(
        #tags$h1("enrichMiR"),
        tags$img(src="enrichMiR_sticker.png"),
        tags$h3("Please wait while the application is initialized..."),
        spin_1()
      )),
      tabItems(
        tabItem(tabName = "tab_intro",
          box(title="enrichMiR: miRNA target enrichment analysis", width=12,
              tags$img(src="enrichMiR_sticker.png", style="float: right;"),
              tags$p(
                "This app will allow you to identify miRNAs whose targets are ",
                "enriched among genesets of interest or a differential ",
                "expression signature, and produce related visualizations.",
                "Although the app was chiefly developed (and benchmarked) for ",
                "miRNAs, some support is also offered to run the same analyses",
                " using RNA-binding proteins."), 
              tags$ol("The app has two main functionalities:", 
                tags$li("performing ", tags$b("target enrichment analysis"),
                        ", either comparing your gene set of interest to a 
                        background set, or using the results of a differential
                        expression analysis (DEA)."),
                tags$li("generating foldchange cumulative distribution (CD) 
                        plots comparing targets and non-targets (requires the
                        results of a differential expression analysis as input).")
              ), tags$br(),
              tags$p("To get started, take a ", 
                     actionLink("helpLink", "quick tour"),
                     " of the app, or browse the help on the upper-right corner."), 
              tags$p("Please report any bug in the ",
                     tags$a(href="https://github.com/ETHZ-INS/enrichMiR", 
                            "github repository"),
                     "."),
              tags$br(),
              tags$p( style="text-align: right;",
                      paste("enrichMiR version",packageVersion("enrichMiR")),
                      "; ", tags$a( href="http://schrattlab.ethz.ch",
                                    "Schratt lab", target="_blank") )
          )
        ),
        
        tabItem(tabName = "tab_species", ############ SPECIES / COLLECTION
          box(title="Select Species and Collection", width=12,
            column(7, selectInput(inputId = "species", "Species", width='98%',
                        choices = c("Human", "Mouse", "Rat", "Fish","Fly","Worm")), tags$br(),
               tags$div(id="collection_input",
                  actionButton(inputId="help_collections", style="float:right;",
                           icon=icon("question-circle"), label=""),
                  selectInput(inputId="collection", width='98%', choices=c(),
                    label="Select a binding sites collection",
                    selected="Targetscan conserved miRNA BS")
                  )),
            column(5, tags$p("Note that some collections (e.g. scanMiR) might 
                             take a moment to load."),
                   withSpinner(textOutput("collection_details")))
          ),
          tags$div(id="exprMirs_box", 
                   box(title="Specify expressed miRNAs (optional)", 
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
                fileInput(inputId="exp_mirna_file", 
                          label="Upload miRNA expression object as '*.csv' file (see below for format)",
                          accept=c("text/csv", ".csv",".tab",".txt",
                                   "text/comma-separated-values,text/plain")),
                div(style = "margin-top: -30px"),
                tags$p("miRNA expression tables should have the following format: miRBase name in the first
                            column, expression-values in the second column"),br(), br(),
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
                    searches will be performed with all miRNAs of the given species"))
        ),
        
        
        tabItem("tab_input", ############ INPUT
          tags$h4("Choose between the two following input options"),
          tags$p("You may either upload the results of a differential expression
                 analysis (DEA), or provide a set of genes of interest against
                 a background set."),
          tabBox(id="input_type", width=12,
            tabPanel(title = "Upload DEA results", value = "dea",
              fluidRow(       
                column(7, fileInput( inputId = "dea_input", 
                  label="Upload DEA object as '*.csv' file (see help)",
                  accept=c("text/csv", ".csv", ".tab", ".txt",
                           "text/comma-separated-values,text/plain")),
                       "Or: ", actionButton("example_dea", "Use example DEA"),
                  tags$hr(), tags$br(),
                  tags$p("Upload a Differential Expression Analysis (DEA) as 
                a table with at least following information: Ensembl ID or Gene 
                Symbol as identifier in the first column, as well as 
                logFC-values and FDR-values."),
                  tags$p("Note that we recommend filtering the DEA to retain
                only the most highly expressed genes (e.g. top 5000).")),
                column(2, actionButton(inputId="help_deaformat",
                          icon=icon("question-circle"), label="View format")),
                column(3, withSpinner(htmlOutput("dea_res")), tags$br(), 
                       tags$hr(), tags$br(),
                     sliderInput(inputId = "dea_sig_th", 
                                 label="Select significance threshold",
                                 min=0.01, max=0.5, value=0.05, step=0.01))
              )
            ),
            tabPanel(title = "Select geneset & background", 
              tags$p("In this mode, your genes of interest are compared against 
                     a background of genes."),
              tags$h3("Genes of interest"),
              tabsetPanel(id="GOI",
                tabPanel(title = "Custom set", value = "GOI_custom",
                         "Paste a list of expressed genes as shown in the window below. A minimum number of two genes is required.",
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
              br(), tags$hr(), tags$h3("Background (required)"), 
              actionButton(inputId="help_background", style="float:right;",
                           icon=icon("question-circle"), label="Help on background"),
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
          )
        ),
        
        
        tabItem(tabName = "tab_enrich", ###### ENRICHMENT ANALYSIS
          column(3, tags$div(id="enrichbtn_outer", 
                             withSpinner(uiOutput("enrichbtn")))), 
          box(width=9, title="Advanced enrichment options", collapsible=TRUE, 
              collapsed=TRUE, sliderInput(inputId="minsize", 
                label="Minium number of annotated targets that is required to
                       consider the miRNA-family for testing",
                          min=1, max=20, value=5, step=1),
              br(),
              box(width="100%", title = "Advanced Options", collapsible=TRUE, collapsed=TRUE, 
                  actionButton(inputId="help_testsadvanced", style="float:right;",
                               icon=icon("question-circle"), label=""),
                  tags$h5(em(
                    "We recommend to only change test settings after reading the",
                    " enrichMiR documentation and the benchmark.")),
                  tags$h5(em("Some tests are always performed by default, namely
                    the 'siteoverlap' test and (except for some annotations and
                    assuming a DEA input) the 'areamir' test.")),
                  uiOutput("extratestinput")
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
                withSpinner(plotlyOutput("bubble_plot")), 
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
                                  choices=ggplot_themes)),
            column(6, downloadLink('bubble_plot_dl', label="Download plot"))
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
          ) )
        ),
        
        
        tabItem(tabName = "tab_cdplot",  ######### CD plot
          box(width=12, title="CD Plot", 
              tags$h3(id="cdplot_na", style="font-color: red;", 
                      icon("exclamation-triangle"), 
                      "Cumulative distribution plots require the use of a DEA input,
                      which you can upload at the input page. You may consult the 
                      tutorial for further info."),
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
                  column(6, textInput("CDplot_xlabel", "x-axis label", 
                                      value="log(foldchange)")),
                  column(6, selectInput("CDplot_theme", "Theme", 
                                        choices=ggplot_themes)),
                  column(6, downloadLink('cd_plot_dl', label="Download plot")),
                  column(6, actionLink('CDplot_dlContent', label="Download plot R data"))
                )
              )
          )
        ),
        tabItem(tabName="tab_benchmark",
          tabBox(width=12,
            tabPanel("Tests description", .testDescription(), 
                     tags$h3("Summary"), 
                     fluidRow(column(12, align="center", tableOutput("testsummary") ))),
            tabPanel("Tests benchmark",
              tags$h3("Benchmark of the different target enrichment tests"),
              tags$p("The different tests were benchmarked on different datasets
                     each involving the transcriptomic characterization of the 
                     knockdown or over-expression of different miRNAs. For each
                     experiment, the signal was additionally scrambled to create 
                     further, more difficult 'pseudo-experiments', which are 
                     averaged in the results below. The benchmark was performed 
                     using TargetScan-predicted sites, and was then used to 
                     guide the choice of default tests in the enrichMiR app."),
              tags$img(src="benchmark1.png"),
              tags$p("Panel A shows the rank of the true miRNA according to 
                     the different tests (lower=better, i.e. a rank of 1 
                     indicates that the true miRNA was correctly identified as
                     the top enriched one). Panel B shows the effective sensitivity 
                     and False Discovery Rate (FDR) of the different tests at a 
                     nominal q-value threshold of 0.05. One can observe that
                     while most tests manage to rank the true hypothesis as 
                     first, most fail to accurately control error."),
              tags$p("In light of these results, the siteoverlap test was 
                     selected as the default for binary signals, and the areamir
                     test for continuous signals. For use with larger annotations
                     (e.g. scanMiR), we however recommend the more conservative
                     lmadd test (see publication for details)."),
              tags$p("Note that restricting the enrichment analysis to the 
                     miRNAs expressed in your system systematically decreases
                     FDR. You can do so in the 'Species and miRNAs' tab, either
                     using a custom list of miRNAs or selecting from available
                     tissues.")
            )
          ), tags$div(style="clear: both;")      
        )
      )
    )
  )
}           




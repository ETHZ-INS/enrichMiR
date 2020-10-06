#' @export
#' 
#' 
#' ######## EnrichMir App Notes
## Species Collection
#0) Which one to take?
#1) We provide Targetscan and KD files for human/mouse/rat
#2) The corresponding miR_family names schould be automatically chosen (at the moment there is an option to choose yourself in the EnrichMir function)
#' 
#' 
## Notes
#1) In general go for the KD-sites. Select that Targetscan can be only done with mice and human?
#2) Think about a colocalization plot
#' 
#' 
#' 
#' 
#' 
enrichMiR.ui <- function(){
  library(shiny)
  library(DT)
  library(shinydashboard)
  library(shinycssloaders)
  library(plotly)
  
  ui <- dashboardPage(
    
    dashboardHeader(title = "EnrichMe", titleWidth = "300px"),
    
    ## Sidebar content
    dashboardSidebar(width = "300px",
                     sidebarMenu(
                       #Either like this or as an own tab
                       selectInput("species", "Species",
                                           choices = c("Human", "Mouse", "Rat","Custom - not yet"), selected = "Human", multiple=FALSE, 
                                            width = '98%'),
                       # menuItem("Select Species",  tabName = "tab_species"),
                       menuItem("Upload Expression Info", 
                                menuSubItem("Upload Background", tabName ="tab_background"),
                                menuSubItem("Optional: Select expressed miRNAs",tabName = "tab_mirnas")
                       ),
                       menuItem("Test Enrichment", 
                                menuSubItem("select collection and genes", "tab_collection"),
                                menuSubItem("enrich", "tab_enrich"),
                                menuSubItem("CD Plot", "tab_cdplot")
                       ),
                       menuItem("Test Colocalization", 
                                menuSubItem("colocalization mode", "tab_co_mode"),
                                menuSubItem("options","tab_co_options"),
                                menuSubItem("colocalize", "tab_colocalize")
                       ),
                       menuItem("About", tabName="tab_about")
                     )
    ),
    ## Body Content
    dashboardBody(
          tabItems(
      #    tabItem(tabName = "tab_species",
      #              tags$h3("Select a Species"), tags$br(),
      #              tabBox(id="species_collection", width=12,
      #                      tabPanel(title="Pre-built", value="prebuilt",
      #                              selectInput("select_speices", "Select Species",
      #                                        choices = c("Human", "Mouse", "Rat"), selected = "Human")),
      #                      tabPanel(title="Upload", value="upload", tags$p("Not yet implemented."))
      #                    ),
      #            box(width=12, withSpinner(verbatimTextOutput("Species_summary",placeholder = TRUE)))
      #            ),
          tabItem(tabName = "tab_background",
                    box(title="Upload the Background", width=12,
                        "Paste a list of expressed genes in 'Ensembl' or 'Gene Symbol' format", br(), br(),
                        textAreaInput("background_genes", 
                                      label="Gene List", 
                                      rows=5,
                                      placeholder="Gene_1\nGene_2\nGene_3", 
                                      resize="vertical"),br(),
                                      footer = "Note: If you want to use the Targetscan miRNA annotations together with rat genes, use the 'Gene Symbol' format")
                  ),
          tabItem(tabName = "tab_mirnas",
                  box(title="Upload expressed microRNAs", width=12,
                      "Paste a list of expressed miRNAs in 'miRBase' format", br(), br(),
                      textAreaInput("expressed_mirnas", 
                                    label="miRNA List", 
                                    rows=5,
                                    placeholder="miRNA_1\nmiRNA_2\nmiRNA_3", 
                                    resize="vertical"),br(),
                      footer = "Note: If no miRNAs are uploaded, enrichment searches will be performed with all microRNAs of the given species")
                  ),
          tabItem(tabName = "tab_collection",
                  tags$h3("Select parameters:"), tags$br(),
                  box(width = 12, title = "Enrichment type",
                    selectInput("collection", "Collection",
                              choices = c("scanMir miRNA BS", "Targetscan miRNA BS", "CISBP RBP motif sites","Custom - not yet"), selected = "scanMir miRNA BS", multiple=FALSE, 
                              width = '98%')
                    ),
                  box(width = 12, title = "Find enrichment in:",
                    tabBox(id="enrichment_type", width=12,
                          tabPanel(title = "Upload DEA", value = "dea",
                                   fileInput(inputId = "dea_input", label = "Upload DEA object as '*.csv' file (see help)",  accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")), br(),
                                   "Upload Differential Expression Analyses (DEAs) as table in the following format: ENSEMBL_ID or Gene Symbol in the first
                                                      column (use same format as for the background), logFC-values in the second column and FDR-values in the third column"
                                   ),
                          tabPanel(title="Custom Genees", value="custom",
                                  textAreaInput("genes_of_interest", 
                                                label="Paste genes of interest", 
                                                rows=5,
                                                placeholder="Gene_A\nGene_B\nGene_C", 
                                                resize="vertical"), br(),
                                  "Use the same annotation format as for the background"
                                  ),
                         tabPanel(title = "GO-Term GeneSet", value = "geneset",
                                  textInput(inputId = "find_go",label = "Find a Go-Term", placeholder = "search with a key word"),
                                  verbatimTextOutput("find_go_result", placeholder = TRUE),
                                  selectizeInput("go-term", "Select Go-Term for enrichment search", choices=c()),))
                    )
                  ),
          tabItem(tabName = "tab_enrich",
                    column(2,actionButton("enrich", "Enrich!", icon = icon("search"))),
                    column(10, tags$h5(textOutput("search for enrichments"))),
                    box(width=12, title="Enrichment Plot", collapsible=TRUE, collapsed=TRUE,
                        withSpinner(plotlyOutput("bubble_plot")),
                        column(6,sliderInput("label.sig.thres","Significance threshold to display labels",min = 0,max = 0.25,value = 0.05,step = 0.01)),
                        column(6,sliderInput("label.enr.thres","Enrichment threshold to display labels",min = 0.5,max = 10,value = 2,step = 0.5)),
                        column(6, numericInput("label_n", "Max number of Labels", value=10, min=1, max=50))
                    ),
                    box(width=12, title="Table", collapsible=TRUE, collapsed = TRUE,
                        withSpinner(DTOutput("hits_table")))
                  ),
          tabItem(tabName = "tab_cdplot",
                    column(6,selectizeInput("mir_fam", "Select miRNA family to display", choices=c(),),
                    actionButton("plotcd", "Plot!", icon = icon("search"))),
                    box(width=12, title="CD Plot", collapsible = TRUE, collapsed = TRUE,
                        withSpinner(plotOutput("cd_plot"))
                        )
                  ),
          tabItem(tabName = "tab_co_mode",
                    tags$h3("Find miRNA colocalizations"), tags$br(),
                    tabBox(id="tabbox_co_mode", width=12,
                         tabPanel(title="In expressed genes", value="expression_coloc",
                                          radioButtons(inputId = "expr_coloc_mode",label = "Choose colocalization type:", choices = c("miRNA - miRNA" = "MM",
                                                                                                                      "miRNA - RBP" = "MR",
                                                                                                                      "miRNA - custom" = "MC"),
                                                      selected = "MM"), br(),
                                          fileInput(inputId = "expr_custominput", label = "Upload custom object as '*.csv' file (see help)",  accept = c(
                                                                                                                                          "text/csv",
                                                                                                                                          "text/comma-separated-values,text/plain",
                                                                                                                                          ".csv"))
                                  ),
                         tabPanel(title = "In subset", value = "subset_coloc",
                                      radioButtons(inputId = "expr_sub_mode",label = "Choose colocalization type:", choices = c("miRNA - miRNA" = "MM",
                                                                                                                           "miRNA - RBP" = "MR",
                                                                                                                           "miRNA - custom" = "MC"),
                                                   selected = "MM"), br(),
                                      fileInput(inputId = "sub_custominput", label = "Upload custom object as '*.csv' file (see help)",  accept = c(
                                        "text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),
                                      textAreaInput("goi_sub_coloc", 
                                                label="Genes of Interest", 
                                                rows=5,
                                                placeholder="Gene_A\nGene_B\nGene_C", 
                                                resize="vertical"),
                                      "Use the same annotation format as for the background", br()
                                  )
                          ),
                  tags$h3("Help text"),
                  column(12,helpText("Include a brief explanation of the two modi. Until now, only Targetscan miRNA sites can be used")),
                  ),
          tabItem(tabName = "tab_co_options",
                  box(title = "Options", width = 12,
                        helpText("Choose one of the following:"),
                        tabBox(id="tabbox_co_dist", width=12,
                             tabPanel(title="Min. Distance", value="min_dist_tab",
                                      sliderInput(inputId = "min_dist",label = "Select a minimum distance between motifs", 
                                                  min = 0, max=250, value = 8,step = 1),
                                      helpText("Minimum number of nucleotides to be located between two motifs that is taken into consideration for the search."), br(),
                                      helpText("By specifying different min / max values, users can search for certain",
                                                "types of functional colocalization. As example, mutual miRNA blocking is",
                                                "suggested to happen within a distance of 0-7Nt. Cooperative repression of",
                                                "two miRNAs is reported to happen within a distance of ca. 8-60Nt.At bigger",
                                                "distances, miRNAs are suggested to act independently on given targets. To search",
                                                "for enrichment of motif pairs on the whole transcript, select 'Minimum Distance", 
                                                "= 0'")
                             ),
                             tabPanel(title="Max. Distance", value="max_dist_tab",
                                      sliderInput(inputId = "max_dist",label = "Select a maximum distance between motifs", 
                                                  min = 0, max=250, value = 7,step = 1),
                                      helpText("Maximum number of nucleotides to be located between two motifs that is taken into consideration for the search."), br(),
                                      helpText("By specifying different min / max values, users can search for certain",
                                               "types of functional colocalization. As example, mutual miRNA blocking is",
                                               "suggested to happen within a distance of 0-7Nt. Cooperative repression of",
                                               "two miRNAs is reported to happen within a distance of ca. 8-60Nt.At bigger",
                                               "distances, miRNAs are suggested to act independently on given targets. To search",
                                               "for enrichment of motif pairs on the whole transcript, select 'Minimum Distance", 
                                               "= 0'")
                             ),
                             tabPanel(title="Distance Range", value="range_tab",
                                      sliderInput(inputId = "range",label = "Select a distance range between motifs", 
                                                  min = 0, max=250, value = c(8,60),step = 1),
                                      helpText("Distance range of nucleotides to be located between two motifs that is taken into consideration for the search."), br(),
                                      helpText("By specifying different min / max values, users can search for certain",
                                               "types of functional colocalization. As example, mutual miRNA blocking is",
                                               "suggested to happen within a distance of 0-7Nt. Cooperative repression of",
                                               "two miRNAs is reported to happen within a distance of ca. 8-60Nt.At bigger",
                                               "distances, miRNAs are suggested to act independently on given targets. To search",
                                               "for enrichment of motif pairs on the whole transcript, select 'Minimum Distance", 
                                               "= 0'")
                             )
                        ), 
                        box(title = "Report option", width = 12,
                          sliderInput(inputId = "min_pairs",label = "Select the minium number of detected pairs displayed in the results:",min = 0,max = 25, step = 1, value = 5)
                        )
                      )
                  ),
          tabItem(tabName = "tab_colocalize",
                    column(2,actionButton("colocalize", "Search for Colocalization!", icon = icon("search"))),
                    column(10, tags$h5(textOutput("Searching for colocalizations. This might take some time!"))),
                    box(width=12, title="Table", collapsible=TRUE, collapsed = TRUE,
                        withSpinner(DTOutput("coloc_table"))),
                    box(width = 12, title = "Some kind of plot",
                        helpText("Either a heatmap and / or a bubble plot or a barplot with the top pairs"))
          )

        )
    )
  )
}           


enrichMir.server <- function(input,output){}          
    
shinyApp(enrichMiR.ui,enrichMir.server)  
    
    
    
    

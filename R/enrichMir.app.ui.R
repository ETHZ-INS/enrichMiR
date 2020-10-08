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
  library(ggplot2)
  
  ui <- dashboardPage(
    
    dashboardHeader(title = "EnrichMe", titleWidth = "300px"),
    
    ## Sidebar content
    dashboardSidebar(width = "300px",
                     sidebarMenu(
                       #Either like this or as an own tab
                       selectInput(inputId = "species", "Species",
                                           choices = c("Human", "Mouse", "Rat","Custom - not yet"), selected = "Human", multiple=FALSE, 
                                            width = '98%'),
                       # menuItem("Select Species",  tabName = "tab_species"),
                       menuItem("Upload Expression Info", 
                                menuSubItem("Upload Background or DEA", tabName ="tab_background"),
                                menuSubItem("Optional: Select expressed miRNAs",tabName = "tab_mirnas")
                       ),
                       menuItem("Select Genes to test", tabName = "tab_set"),
                       menuItem("EnrichMir", 
                                menuSubItem("Enrichment Options", tabName = "tab_options"),
                                menuSubItem("enrich", "tab_enrich"),
                                menuSubItem("CD Plot", "tab_cdplot")
                       ),
                       menuItem("CoMir", 
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
                    box(title="Upload Background or DEA", width=12,
                      tabBox(id="upload_background", width=12,
                             tabPanel(title = "Upload Background", value = "upload",
                                      "Paste a list of expressed genes in 'Ensembl' or 'Gene Symbol' format", br(), br(),
                                      textAreaInput(inputId = "background_genes", 
                                                    label="Gene List", 
                                                    rows=5,
                                                    placeholder="Gene_1\nGene_2\nGene_3", 
                                                    resize="vertical"),br(),
                                      footer = "Note: If you want to use the Targetscan miRNA annotations together with rat genes, use the 'Gene Symbol' format"),
                             tabPanel(title = "Upload DEA", value = "dea",
                                      fileInput(inputId = "dea_input", label = "Upload DEA object as '*.csv' file (see help)",  accept = c(
                                        "text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")),
                                      checkboxInput("header", "Header", TRUE), br(),
                                      "Upload Differential Expression Analyses (DEAs) as table in the following format: ENSEMBL_ID or Gene Symbol in the first
                                                      column (use same format as for the background), logFC-values in the second column and FDR-values in the third column"
                             ))
                    )
                  ),
          tabItem(tabName = "tab_mirnas",
                  box(title="Upload expressed microRNAs", width=12,
                      "Paste a list of expressed miRNAs in 'miRBase' format", br(), br(),
                      textAreaInput(inputId = "expressed_mirnas", 
                                    label="miRNA List", 
                                    rows=5,
                                    placeholder="miRNA_1\nmiRNA_2\nmiRNA_3", 
                                    resize="vertical"),br(),
                      footer = "Note: If no miRNAs are uploaded, enrichment searches will be performed with all microRNAs of the given species")
                  ),
          tabItem(tabName = "tab_set",
                  tags$style(HTML("
                      .tabbable > .nav > li[class=active]    > a {background-color: dimgrey; color:white}
                      .tabbable > .nav > li > a                  {background-color: gray75;  color: black}
                                      ")),
                  box(width = 12, title = "Choose Gene Set:",
                      tabsetPanel(id="choose_gene_set", 
                                  tabPanel(title = "Custom Geen Sets", value = "custom_gene_sets",
                                           br(),
                                           tags$h4("Select parameters:"),
                                           tabBox(id="cust_enrichment_type", width=12,
                                                  tabPanel(title="Genes of Interest", value="custom_genes",
                                                           textAreaInput(inputId = "genes_of_interest", 
                                                                         label="Paste genes of interest", 
                                                                         rows=5,
                                                                         placeholder="Gene_A\nGene_B\nGene_C", 
                                                                         resize="vertical"), br(),
                                                           "Use the same annotation format as for the background"
                                                  ),
                                                  tabPanel(title = "GO-Term GeneSet", value = "go_genes",
                                                           fluidRow(
                                                           box(title = "Find a Go-Term",width = 12, collapsible = TRUE,collapsed = TRUE,
                                                           textInput(inputId = "find_go",label = "Type in keyword", placeholder = "search with a key word"),
                                                           column(6,radioButtons(inputId = "ontology",label = "Ontology:", choices = c("Biological Process" = "BP",
                                                                                                                                       "Cellular Component" = "CC",
                                                                                                                                       "Molecular Function" = "MF"),
                                                                                 selected = "CC")), br(),
                                                           column(6,actionButton(inputId = "find_go_button", "Find GO!", icon = icon("search"))),
                                                           br(),
                                                           column(12,"The first 20 results of the search are displayed"),
                                                           column(12,verbatimTextOutput("find_go_result", placeholder = TRUE)),
                                                           ),
                                                           column(12,radioButtons(inputId = "go_genes_format",label = "Background Gene Format:", choices = c("Ensembl" = "Ens",
                                                                                                                                                      "Gene Symbol" = "GS"),
                                                                                  selected = "Ens")),
                                                           column(12,selectizeInput(inputId = "go_term", "Select Go-Term for enrichment search", choices=c()))
                                                           )
                                                  ))
                                           
                                  ),
                                  tabPanel(title = "DEA Sets", value = "dea_sets",
                                           br(),
                                           tags$h4("Set threshhold for binary enrichment tests:"), br(),
                                           sliderInput(inputId = "binary_fdr", label = "Select a FDR cut-off",
                                                       min = 0.05,max = 0.5,value = 0.05, step = 0.05,width = "80%")
                                  )
                      )) 
                  ),
          tabItem(tabName = "tab_options",
                  tags$h3("Select enrichment options:"), tags$br(),
                  box(width = 12, title = "Collections",
                    selectInput("collection", label = "",
                              choices = c("scanMir miRNA BS", "Targetscan miRNA BS", "CISBP RBP motif sites","Custom - not yet"), selected = "scanMir miRNA BS", multiple=FALSE, 
                              width = '98%'),
                    "Until now only Targetscan"
                    ),
                  tabBox(id="test_type", width=12,
                         tabPanel(title = "Binary Test", value = "binary",
                                  sliderInput(inputId = "minsize", label = "Select the minium number of targets to be considered for testing",
                                              min = 1,max = 20,value = 5, step = 1),
                                  br(),
                                  "Brief description of binary test"),
                         tabPanel(title = "Continuous Test",value = "continous",
                                  radioButtons(inputId = "cont_test_buttons", label = "Choose to perform a continuous enrichment test based on:",
                                               choices = c("Site Scores (only miRNAs)" = "SC",
                                                           "logFC" = "LFC"), 
                                               selected = "SC"),
                                  br(),
                                  "In continuous testing mode, all genes get used for enrichment analysis. Brief explanation of continuous test")
                         
                        )
                  ),
          tabItem(tabName = "tab_enrich",
                    column(2,actionButton(inputId = "enrich", "Enrich!", icon = icon("search"))),
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
                    column(12,selectInput(inputId = "mir_fam", "Select miRNA family to display", choices=c())),
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
                                        ".csv"))
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
                    column(2,actionButton(inputId = "colocalize", "Search for Colocalization!", icon = icon("search"))),
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


enrichMiR.server <- function(input,output){}          
    
shinyApp(enrichMiR.ui,enrichMiR.server)  
    
    
    
    

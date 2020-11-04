

#' @import shiny DT shinydashboard shinycssloaders
#' @export
enrichMiR.ui <- function(){
  library(shiny)
  library(shinyjqui)
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
                       menuItem("Species and miRNAs", tabName = "tab_species"),
                       menuItem("Input genes/DEA", tabName = "tab_input"),
                       menuItem("Enrichment analysis",
                                menuSubItem("Enrichment Options", tabName = "tab_enrich_op"),
                                menuSubItem("Enrich!", tabName = "tab_enrich"),
                                menuSubItem("CD Plot", tabName = "tab_cdplot")
                       ),
                       menuItem("Co-localization", 
                                menuSubItem("colocalization mode", tabName = "tab_co_mode"),
                                menuSubItem("options",tabName = "tab_co_options"),
                                menuSubItem("colocalize", tabName = "tab_colocalize")
                       ),
                       menuItem("About", tabName="tab_about")
                     )
    ),
    
    ## Body Content
    dashboardBody(
      tags$head(tags$style(HTML("
#sel_test_div label {
  display: table-cell;
  text-align: center;
  vertical-align: middle;
}
#sel_test_div .form-group {
  display: table-row;
}
"))),
      tabItems(
        tabItem(tabName = "tab_species",
                    box(title = "Select Species and Collection", width=12,
                    selectInput(inputId = "species", "Species", width = '98%', multiple=FALSE, 
                                choices = c("Human", "Mouse", "Rat","Custom - not yet"), 
                                selected = "Human"), tags$br(),
                    selectInput(inputId="collection", label = "Select a binding sites collection", width = '98%',
                                choices = c("scanMir miRNA BS", "Targetscan conserved miRNA BS","Targetscan all miRNA BS", "CISBP RBP motif sites","miRTarBase", "Custom - not yet"),
                                selected = "Targetscan conserved miRNA BS", multiple=FALSE)
                    ),
                box(title = "Expressed miRNAs", collapsible=TRUE, collapsed=TRUE, width=12,
                    tabBox(id="expressed_mirna_box", width=12,
                    tabPanel(title = "Custom Set",
                             "Paste a list of expressed miRNAs in 'miRBase' format", br(), br(),
                              textAreaInput(inputId = "exp_mirna_list", 
                                  label="miRNA List", 
                                  rows=5,
                                  placeholder="miRNA_1\nmiRNA_2\nmiRNA_3", 
                                  resize="vertical")),
                    tabPanel(title = "Upload miRNA expression table",
                             "Upload miRNA Expression Table in the following format: miRBase name in the first
                                          column, expression-values in the second column", br(),
                             fileInput(inputId = "exp_mirna_file", label = "Upload miRNA expression object as '*.csv' file (see help)",  accept = c(
                                  "text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
                                checkboxInput("header_mir", "Header", TRUE),
                                br(),
                                sliderInput(inputId = "mir_cut_off", label = "Select miRNA expression cut-off:",
                                            min = 10, max = 100, post  = " %", value = 50 ),
                                "Keep only the top ..% expressed miRNAs"
                             )),
                    br(),
                    footer = "Note: If no miRNAs are uploaded, enrichment searches will be performed with all microRNAs of the given species")
                # eventually enable this as upload or selection from pre-loaded tissues?
        ),
        tabItem("tab_input",
                tabBox(id="input_type", width=12,
                       tabPanel(title = "Select geneset & background", 
                                tags$p("In this mode, your genes of interest are compared against a background of genes."),
                                tags$h3("Gene Format"),
                                radioButtons(inputId = "genes_format", label = "Select:", 
                                                 choices = c("Ensembl" = "Ens",
                                                             "Gene Symbol" = "GS"),
                                                 selected = "Ens"),
                                br(),tags$hr(),
                                tags$h3("Genes of interest"),
                                tabsetPanel(id="GOI",
                                            tabPanel(title = "Custom set", value = "GOI_custom",
                                                     "Paste a list of expressed genes in the selected format", br(), br(),
                                                     textAreaInput(inputId = "genes_of_interest", 
                                                                   label="Gene List", 
                                                                   rows=5,
                                                                   value = NULL,
                                                                   placeholder="Gene_1\nGene_2\nGene_3", 
                                                                   resize="vertical"),
                                            ),
                                            tabPanel(title = "From Gene Ontology", value="GOI_GO",
                                                     # selectizeInput(inputID="go_type", "Select GO collection", 
                                                     #                choices=c("Cellular component"="CC",
                                                     #                          "Biological process"="BP",
                                                     #                          "Molecular function"="MF")),
                                                     selectizeInput(inputId="go_term", "Select GO-Term", choices=c()),
                                                     textOutput("GOI_nb"))
                                ),br(),tags$hr(),tags$h3("Background"),
                                textAreaInput(inputId = "background_genes", 
                                              label="Gene List", 
                                              rows=5,
                                              value = NULL,
                                              placeholder="Gene_1\nGene_2\nGene_3", 
                                              resize="vertical"),br(),
                                footer = "Note: If you want to use the Targetscan miRNA annotations together with rat genes, use the 'Gene Symbol' format"
                       ),
                       tabPanel(title = "Upload DEA results", value = "dea",
                                fileInput(inputId = "dea_input", label = "Upload DEA object as '*.csv' file (see help)",  accept = c(
                                  "text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
                                checkboxInput("header", "Header", TRUE), br(),
                                sliderInput(inputId = "dea_sig_th", label = "Select significance thresshold",
                                            min = 0.01, max = 0.5, value = 0.05, step = 0.01),
                                "Upload Differential Expression Analyses (DEAs) as table with at least following information: Provide ENSEMBL_ID or Gene Symbol as identifier
                                in the first column (same format as for the background), as well as logFC-values and FDR-values"
                       )
                )
        ),
        tabItem(tabName = "tab_enrich_op",
                box(title = "Select enrichment options:", collapsible=TRUE, collapsed=FALSE, width=12,
                    sliderInput(inputId = "minsize", label = "Minium number of annotated targets that is required to consider the miRNA-family 
                                for testing",
                                min = 1,max = 20,value = 5, step = 1),
                    br(),
                    box(title = "Advanced Options", collapsible=TRUE, collapsed=TRUE, width=12,
                        tags$h5(em("We recommend to only change test settings after reading the EnrichMir manual and the Benchmark descriptions in the paper")),
                        tags$h5(em("By default, we always perform the 'siteoverlap' and the 'areamir' tests")),
                        checkboxGroupInput(inputId = "tests2run", label = "Select additional tests to run", choices = list(
                                                                                                                    "overlap" = "overlap",
                                                                                                                    "weight.overlap" ="woverlap",
                                                                                                                    "plMod" = "plMod",
                                                                                                                    "modscore" = "modscore",
                                                                                                                    "ks" = "ks",
                                                                                                                    "mw" = "mw"), selected = NULL, inline = FALSE,width = NULL)
                        )
                    )
                ),
        tabItem(tabName = "tab_enrich",
                column(6,actionButton(inputId = "enrich", "Enrich!", icon = icon("search"))),
                column(4, id="sel_test_div", selectInput("view_test", "Test:", choices=c(), width="450px")),
                column(2, checkboxInput("view_all", "advanced")),
                
                # box(width=12, title="Select test to view", collapsible=TRUE, collapsed=TRUE,
                #     tabBox(id="test_type", width=12, 
                #            tabPanel(title = "Binary Test", value = "binary",
                #                     selectInput("view_binary", "Interested in up- or downregulated genes:", choices = c("")),
                #                     # radioButtons(inputId = "up_down", label = "Interested in up- or downregulated genes:",
                #                     #              choices = c("Up" = ".up",
                #                     #                          "Down" = ".down"),
                #                     #              selected = ".down"),
                #                     "A hypergeometric test on the number of binding sites will be performed to calculate significance"),
                #            tabPanel(title = "Continuous Test",value = "continous",
                #                     tags$h4(em("Only with DEA-Input")),br(),
                #                     br(),
                #                     "In continuous testing mode, all genes get used for enrichment analysis.Here, we employ an analytic 
                #                     rank-based enrichment analysis using a conversion of the scores as weights."),
                #            tabPanel(title = "Advanced - Select test by name",value = "advanced",
                #                     tags$h4(em("In case you've selected additional tests for the Enrichment Analyses,
                #                                you can select them here. Again, please refer to the manual as 
                #                                well as the benchmark results for individual test infos!")),br(),
                #                     selectInput("view_test", "View test", choices = c(""))
                #            )
                #     ),
                #     tags$h5("Name of chosen test:"),
                #     verbatimTextOutput(outputId = "test_info"),
                #     ),
                box(width=12, title="Enrichment Plot", collapsible=TRUE, collapsed=TRUE,
                    withSpinner(jqui_resizable(plotlyOutput("bubble_plot"))),
                    column(6,sliderInput(inputId = "label.sig.thres","Significance threshold to display labels",min = 0,max = 0.25,value = 0.05,step = 0.01)),
                    column(6,sliderInput(inputId = "label.enr.thres","Enrichment threshold to display labels",min = 0.5,max = 10,value = 1,step = 0.5)),
                    column(6, numericInput(inputId = "label_n", "Max number of Labels", value=10, min=1, max=50)),
                    column(6, radioButtons(inputId = "sig.field", label = "Display on y-axis:",
                                           choices = c("p.value" = "pvalue",
                                                       "FDR" = "FDR"), selected = "FDR"))
                ),
                box(width=12, title="Table", collapsible=TRUE, collapsed = TRUE,
                    checkboxGroupInput(inputId = "columns2show", label = "Select add. columns to be shown:", choices = list(
                      "miRNA names" = "members",
                      "predicted target genes" ="genes"), selected = NULL, inline = TRUE),
                    withSpinner(DTOutput("hits_table")),
                    downloadLink('dl_hits', label = "Download all"))
        ),
        tabItem(tabName = "tab_cdplot",
                box(width=12, title="CD Plot", 
                    "(CD plots require a DEA input.)",
                    br(), br(),
                    column(6,selectizeInput(inputId = "mir_fam", "Select miRNA family to display", choices=c())),
                    column(6,selectInput(inputId = "CD_type", "Split by", choices=c("sites","score"))),
                    withSpinner(jqui_resizable(plotOutput("cd_plot",width = '100%', height = '400px'))),
                    br(),br(),br(),br(),br(),
                    column(12,sliderInput(inputId = "CDplot_xaxis","logFC to display on x.axis",min = 0.5,max = 5,value = 2,step = 0.5))
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




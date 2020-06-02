enrichMiR.ui <- function(){
  library(shiny)
  library(shinydashboard)
  
  ## Should we give the minDist option in the TopTarget Search?
  ## Name equal inputs in different Tabs equally or differently?
  ## Do we need the Update Species Button
  ## Include a filter KD slider? Either in Advanced Search options, or in Output Boxes?
  
  
  
  ui <- dashboardPage(
    
    dashboardHeader(title = "findMirSite", titleWidth = "300px"),
    
    ## Sidebar content
    dashboardSidebar(width = "300px",
                     sidebarMenu(
                       #Basic Infos
                       selectInput(inputId = "Genome", label = "Select an Ensembl Genome", choices = list("Mouse Mus_Musculus (mm10)", 
                                                                                                          "Rat Rattus_Norvegicus (rnor6)",
                                                                                                          "Human Homo_Sapiens (hg38)"),
                                   selected = "Mouse Mus_Musculus (mm10)"),
                       selectInput(inputId = "miRNA_Species", label = "Select an miRNA Species", choices = list( "Mouse", 
                                                                                                                 "Rat",
                                                                                                                 "Human"),
                                   selected = "Mouse"),
                       menuItem("Custom Search", tabName = "CustomSearch", icon = icon("user-edit")),
                       menuItem("Gene/Transcript Search", tabName = "GeneSearch", icon = icon("book-reader")),
                       menuItem("miRNA Top Targets", tabName = "miRNATargets", icon = icon("book")),
                       menuItem("Plot miRNA affinity", tabName = "AffinityPlot", icon = icon("chart-bar"))
                       #actionButton(inputId = "UpdateSpecies", label = "Update Species")
                     )
    ),
    ## Body Content
    dashboardBody(
      tabItems(
        # Custom Search tab content
        tabItem(tabName = "CustomSearch",
                fluidRow(
                  column(width = 4,
                         
                         ##Input
                         #Paste sequence
                         box(title = "Paste Sequence", width = NULL,
                             textAreaInput(inputId = "PasteSequence_Custom-Search",label = "RNA or DNA format", placeholder = "Paste in here a sequence you want to search for binding sites",
                                           height = "200px"),
                             checkboxInput(inputId = "is-circular_Custom-Search", label =  "is circular")
                         ),
                         #miRNA selection
                         tabBox(title = "Search Options", side = "right", selected = "miRNA", width = NULL,
                                tabPanel(title = "Advanced",
                                         sliderInput(inputId = "minDist_Custom-Search", label = "Minium distance in NT between two binding sites of the same miRNA", 
                                                     min = 0, max = 7, value = 0, step = 1),
                                         checkboxInput(inputId = "KeepMatchSeq_Custom-Search", label =  "Keep matched sequence")
                                ),
                                tabPanel(title = "miRNA",
                                         radioButtons(inputId = "miRNA_select_Custom-Search",label = "Search for:", choiceNames = list("A specific miRNA",
                                                                                                                                       "All miRNAs",
                                                                                                                                       "Conserved miRNAs in mammals"),
                                                      choiceValues = list("A specific miRNA_C",
                                                                          "All miRNAs_C",
                                                                          "Conserved miRNAs in mammals_C")
                                         ),
                                         textInput(inputId = "type-in-miRNA_Custom-Search", label = "miRNA name in mirbase format",
                                                   placeholder = "e.g. mmu-miR-134-5p")
                                )),
                         #Action Button
                         box(width = NULL,
                             actionButton(inputId = "SearchButton_Custom", label = "Search")
                         )),
                  column(width = 7,
                         ##Output
                         box(title = "ggbio_ManhattanPlot_Output", width = NULL,
                             plotOutput(outputId = "Manhattan_Custom-Search")
                         ),
                         box(title = "DataTable-Output", width = NULL,
                             dataTableOutput(outputId = "DataTable_Custom-Search"),
                             downloadButton(outputId = "DownloadDT_Custom-Search", label = "Download")
                         ))
                )),
        
        # Gene Search tab content
        tabItem(tabName = "GeneSearch",
                fluidRow(
                  column(width = 4,
                         ##Input
                         #Gene Selection
                         box(title = "Gene Info", width = NULL,
                             selectInput(inputId = "Gene-Selection_Gene-Search", choices = list("Gene Name" = "Gene_Name", 
                                                                                                "Ensembl Gene ID" = "EnsGene",
                                                                                                "Ensembl Transcript ID" = "EnsTrans"),
                                         label = "Choose the gene or transcript 3'UTR to search for binding sites", selected = "Gene_Name"),
                             textInput(inputId = "type-in-Gene_Gene-Search", label = NULL, placeholder = "Limk1")
                         ),
                         #miRNA selection
                         tabBox(title = "Search Options", side = "right", selected = "miRNA", width = NULL,
                                tabPanel(title = "Advanced",       
                                         sliderInput(inputId = "minDist_Gene-Search", label = "Minium distance in NT between two binding sites of the same miRNA", 
                                                     min = 0, max = 7, value = 0, step = 1),
                                         sliderInput(inputId = "Shadow_Gene-Search", label = "Ribosomal shadow at the beginning of the 3'UTR, recommended is 15 NT", 
                                                     min = 0, max = 20, value = 15, step = 1),
                                         checkboxInput(inputId = "KeepMatchSeq_Gene-Search", label =  "Keep matched sequence")
                                ),
                                tabPanel(title = "miRNA",
                                         radioButtons(inputId = "miRNA_select_Gene-Search",label = "Search for:", choiceNames = list("A specific miRNA",
                                                                                                                                     "All miRNAs",
                                                                                                                                     "Conserved miRNAs in mammals"),
                                                      choiceValues = list("A specific miRNA_G",
                                                                          "All miRNAs_G",
                                                                          "Conserved miRNAs in mammals_G")
                                         ),
                                         textInput(inputId = "type-in-miRNA_Gene-Search", label = "miRNA name in mirbase format",
                                                   placeholder = "mmu-miR-134-5p")
                                )),
                         box(width = NULL,                  
                             #Action Button
                             actionButton(inputId = "SearchButton_Custom", label = "Search")
                         )),
                  column(width = 7,
                         ##Output
                         box(title = "ggbio_ManhattanPlot_Output", width = NULL,
                             plotOutput(outputId = "Manhattan_Gene-Search")
                         ),
                         box(title = "DataTable-Output", width = NULL,
                             dataTableOutput(outputId = "DataTable_Gene-Search"),
                             downloadButton(outputId = "DownloadDT_Gene-Search", label = "Download")
                         )
                  ))),
        
        
        # Top Target tab content
        tabItem(tabName = "miRNATargets",
                fluidRow(
                  column(width = 4,
                         ##Input
                         #miRNA selection
                         tabBox(title = "Search Options", side = "right", selected = "miRNA", width = NULL,
                                tabPanel(title = "Advanced",  
                                         sliderInput(inputId = "minDist_Top-Target", label = "Minium distance in NT between two binding sites of the same miRNA", 
                                                     min = 0, max = 7, value = 0, step = 1),
                                         sliderInput(inputId = "Shadow_Top-Target", label = "Ribosomal Shadow at the beginning of the 3'UTR, recommended is 15 NT", 
                                                     min = 0, max = 20, value = 15, step = 1)
                                ),
                                tabPanel(title = "miRNA",        
                                         textInput(inputId = "type_in-miRNA_Top-Target", label = "miRNA name in mirbase format",
                                                   placeholder = "mmu-miR-134-5p")
                                )),
                         box(title = "Future Options", width = NULL,
                             "it would be cool to for example subset for", br(), "GO-Terms in the future"
                         ),
                         #Action Button
                         box(width = NULL,
                             actionButton(inputId = "SearchButton_Custom", label = "Search"))
                  ),
                  
                  column(width = 7,
                         ##Output
                         box(title = "DataTable-Output", width = NULL,
                             dataTableOutput(outputId = "DataTable_Top-Target"),
                             downloadButton(outputId = "DownloadDT_Top-Target", label = "Download")
                         )
                  )
                )),
        
        
        # Affinity Plot tab content
        tabItem(tabName = "AffinityPlot",
                fluidRow(
                  
                  ##Input
                  #miRNA selection
                  box(width = 4,
                      textInput(inputId = "type_in-miRNA_Affinity", label = "miRNA name in mirbase format",
                                placeholder = "mmu-miR-134-5p"),
                      
                      #Action Button
                      actionButton(inputId = "PlotButton_Affinity", label = "Plot")
                  ),
                  
                  ##Output
                  box(title = "KdPlot_Output", width = 7,
                      plotOutput(outputId = "KDPlot_Affinity"),
                      downloadButton(outputId = "DownloadKDPlot_Affinity", label = "Download")
                  )
                ))
      )
    ))
}
#' enrichMiR.server
#'
#' @return A shiny server function
#' @export
enrichMiR.server <- function(){
  library(DT)
  library(ggplot2)
  
  dtwrapper <- function(d, pageLength=25){
    datatable( d, filter="top", class="compact", extensions=c("Buttons","ColReorder"),
               options=list(pageLength=pageLength, dom = "fltBip", rownames=FALSE,
                            colReorder=TRUE, 
                            buttons=c('copy', 'csv', 'excel', 'csvHtml5') ) )
  }
  
  
  function(input, output, session){
    
    ##############################
    ## initialize inputs
    
    # Here has to come the Go-Term Inputs
    # And the miRNA family input
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
}

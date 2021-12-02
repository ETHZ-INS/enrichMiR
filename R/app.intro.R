.getAppIntro <- function(){
  data.frame(
    element=c(NA, "#menu_species", "#collection_input", "#exprMirs_box", 
              "#expressed_mirna_outer", "#menu_input", "#input_type",
              "#example_dea", "#menu_input", "#menu_enrich", "#enrichbtn_outer",
              "#resultsbox", "#menu_cdplot", "#mir_fam", "#cdplot_outer", NA),
    intro=c(
      "This introduction will walk you through the usage of the enrichMiR app.
      <br/><br/>
      You can use the Next/Previous buttons to navigate the tutorial, and leave
      it at any time by clicking outside it.",
      "The first step is to make sure that you're working with the right species
      and annotation. To set this, you open the 'Species and miRNAs' tab.",
      "Once you've selected a species, you will be able to see and select one
       of its available target annotations (see the help button for information
       about the annotations). For the sake of this example, leave it at the 
      starting default, i.e. human 'targetScan conserved miRNA BS'.",
      "This is optional, but you can also specify expressed miRNAs by expanding
       this box (click the \"+\" on the right).<br/>This will be used to 
       restrict the miRNAs considered for enrichment analysis. In addition, you
       will be able to visualize given expression values in the results.",
      "There are three ways to provide this information:<br/>
       1) you can specify a list of expressed miRNAs in 'Custom set';<br/>
       2) you can upload a table of miRNA expression levels;<br/>
       3) you can pick pre-compiled miRNA expression profiles for your 
       tissue/celltype of interest.<br/><br/>Not that you do not necessarily
       need to specify expressed miRNA, but this is likely to give you better
       results.",
      "Once this is done, the next step is to provide the signal in which 
       enrichment should be looked for (e.g. your genes of interest, or 
       differential expression signature). To do this, click on the Input tab.",
      "There are two main types of input:<br/>
      1) the results of a differential expression analysis (DEA), which you can
      upload using the 'Upload DEA results' tab, or<br/>
      2) a set of genes to be tested for over-representation (the set can be 
      specified either manually or by selecting genes belonging to a given GO 
      category)<br/>
      When specifying a set of genes, these will be tested for 
      over-representation among targets in comparison to a background set of 
      genes (also known as 'universe'). This is critical to over-representation
      analysis! See the corresponding help button for more information about the
      choice of a background.<br/><br/>
      For the purpose of this example, select the 'Upload DEA results' tab.",
      "For this tutorial, we'll use an in-built example DEA, which you can load
      by clicking the 'Use example DEA' button.<br/><br/>Once you've done so,
      you will see on the right that it successfully read the DEA file.",
      "If you've provided a valid input, the input menu on the left will show
      a blue badge instead of a red one.",
      "We're now ready to launch the enrichment analysis!<br/><br/>First go to
      the 'Enrichment analysis' tab",
      "Then click the 'Run enrichMir!' button. The analysis will take a moment.",
      "Once it's over, an enrichment plot will appear, displaying the 
      results. You can hover on the points for more information on the top hits.
      (For more information on the enrichment plot and its interpretation, click
      the help button above it).<br/><br/>
      You can also visualize (and download) the results as a table in the 
      'Results Table' tab.",
      "Independently of an enrichment analysis, when your input is a 
      differential expression analysis (DEA) results you also have the 
      possibility to create cumulative distribution (CD) plots for predicted 
      targets of any miRNA. To try this out, open the 'CD plot' tab.",
      "By default, a CD plot is created for the first miRNA in the alphabetical
      list. To select your miRNA of interest, click in the selector, erase the 
      current miRNA (using the backspace button), and start typing, for example
      '133a'. The miRNA compatible with your selection will show up, and you can
      select the one you are interested in. Select a miR-133a-3p, which is a 
      member of the family whose targets showed the strongest enrichment in 
      this DEA.",
      "After a moment, the cumulative distribution plot will appear, showing
      the foldchange distributions of different types of targets and 
      non-targets. In the presence of a real miRNA effect, we expect to see a 
      gradual dose-response pattern (targets with stronger binding sites also 
      have stronger log-foldchanges), as is the case in this example.",
      "This concludes our tutorial. While there are additionnal options here 
      and there, these are the app's main functionalities."
    )
  )
}

.getHelpModal <- function(topic){
  switch(topic,
         collections=modalDialog(title="Binding sites collections", tags$p(
"Each collection represents a set of binding site (validated or predicted) for
(a set of) miRNAs or families. Not all collections are available in all species.
Here is a brief description of each type of collection:"),
           tags$ul(
             tags$li(tags$b("TargetScan: "),
               "miRNA binding site predictions from ",
               tags$a( href="http://www.targetscan.org","TargetScan",
                       target="_blank"),
               " (version 8). These are available in two flavors: using all 
               sites, or only sites conserved across vertebrates. In both cases,
               only one UTR per gene (typically the longest) is considered."),
             tags$li(tags$b("ScanMiR: "),
               "miRNA binding site predictions using ",
               tags$a(href="https://www.bioconductor.org/packages/devel/bioc/html/scanMiR.html",
                      "scanMiR", target="_blank"),
               "a transcript-level affinity prediction method."),
             tags$li(tags$b("oRNAment: "),
               "RNA binding protein (RBPs) bindings from the ",
               tags$a(href="http://rnabiology.ircm.qc.ca", "oRNAment", 
                      target="_blank"), " database.", tags$br(),
               "Note that enrichMiR was developed for miRNA bindings, and its 
               usage for RBPs is experimental.")
            )
         ),
         enrichplot=modalDialog(title="Enrichment plot", tags$p(
           "The enrichment plot shows the significance of the enrichment on the
           y-axis, and the magnitude of the enrichment on the x-axis. For basic
           overlap/over-representation tests, the enrichment is the fold 
           increase over the overlap expected by chance. For continuous tests
           (based on the results of a differential expression analysis), the 
           enrichment represents a magnitude of association with the 
           log-foldchanges of the putative targets.<br/>The size of the points 
           is proportional to the number of targets with binding site for the 
           given miRNA/family. If miRNA expression data was given, the points 
           are colored by this expression."),
           tags$p("Hovering on the points/families will provide extra 
                  information about the family members and the test statistics.")
         ),
         cdplot=modalDialog(title="Cumulative distribution (CD) plots", tags$p(
           "Cumulative distribution plots are used to visualize differences in
           log-foldchange distributions between groups of genes/transcripts.
           Briefly, the genes in each sets are ranked according to their 
           foldchange (the x-axis), while the y-axis shows the proportion of 
           genes in the set that have a foldchange below or equal to that given
           on x. In this way, small shifts in foldchanges can easily be 
           visualized."), tags$p(
           "There are different ways of splitting the genes into groups, which
           can be selected using the 'Split by' dropdown list. We recommend 
           splitting by best binding type if this information is available in 
           the collection you are using, or otherwise splitting by predicted
           repression score. This is the default (Automatic) procedure."),
           tags$p("In the presence of a real miRNA effect, we expect to see a 
            gradual dose-response pattern, meaning that targets with stronger 
            binding sites also have stronger log-foldchanges.")
         ),
         modalDialog(title=topic,
                     "No help currently available for this topic.")
  )
}
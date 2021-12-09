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
               "a transcript-level affinity prediction method. This target 
               annotation is compatible with both gene and transcript IDs."),
             tags$li(tags$b("miRTarBase: "),
                     "Experimentally-validated miRNA targets from the ",
                     tags$a(href="https://mirtarbase.cuhk.edu.cn/~miRTarBase/",
                            "miRTarBase", target="_blank"),
                     "database. While these targets are more trustworthy, a 
                     sufficient number of targets is available for only a 
                     proportion of the miRNAs, and some targets might have been
                     missed due to the cellular context of the experiments. 
                     Because this collection does not include individual sites,
                     only the simplest over-representation test is available."),
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
           repression score (if available). This is the default (Automatic) 
           procedure."),
           tags$p(renderPlot({
             CDplot(list("no site"=rnorm(600,mean=0,sd=0.8), 
                         "7mer"=c(rnorm(200,0),rnorm(200, mean=-0.5)), 
                         "8mer"=c(rnorm(100,0),rnorm(200, mean=-1))), size=1.3) + 
               scale_color_manual(values=enrichMiR:::.siteTypeColors()[
                 c("no site", "7mer", "8mer")]) + xlab("logFC") + xlim(-3,2)
           })),
           tags$p("In the presence of a real miRNA effect, we expect to see a 
            gradual dose-response pattern, as shown above, meaning that 
            targets with stronger binding sites also have stronger 
            log-foldchanges (i.e. lower logFC, in the case of an increased 
            miRNA activity, as in this example)."),
           tags$p("If your CD plot shows very angular lines, rather than 
            lines that approach curves, this means that your sets do not contain
            enough elements. An example of this is given below:"),
           tags$p(renderPlot({
             CDplot(list("no site"=rnorm(600,mean=0,sd=0.8), 
                         "7mer"=c(rnorm(8,0),rnorm(15, mean=-0.5)), 
                         "8mer"=c(-2.8, -0.5, 0.5)), size=1.3) + 
               scale_color_manual(values=enrichMiR:::.siteTypeColors()[
                 c("no site", "7mer", "8mer")]) + xlab("logFC") + xlim(-3,2)
           })),
           tags$p("In such cases, if you are using score intervals, try reducing
            the number of sets in the plot options. If you are using sites or 
            site types, this means that the miRNA of interest does not have 
            enough sites of this type in your dataset. This can especially 
            happen when using TargetScan ", tags$em("conserved"), "sites, in 
            which case you can try using all predicted sites instead.")
         ),
         tests=modalDialog(title="Enrichment tests", tags$p(
            "enrichMiR implements different statistical tests for target 
            enrichment. Some tests depend on continous inputs (and hence require
            a DEA analysis to be provided), while others are binary, based on 
            the membership in a gene set. When the input is a DEA, each binary
            tests will be performed on both significantly upregulated and 
            downregulated genes."
           ), tags$p(
             "Some tests gave poor performances in the benchmark, and should 
             therefore not be used. By default, the 'siteoverlap' and 'areamir' 
             tests are enabled -- these are those that gave the best performance."
           ), tags$p(
             "For more information on the tests, please see the documentation."
           ), tags$p(
             "Finally, note that the tests were developed and benchmarked for
             application with predicted miRNA targets, not RNA binding proteins
             (RBPs), which have more degenerate binding patterns. Although we 
             offer the possibility to perform RBP target enrichment analyses, 
             the tests were not benchmarked for that purpose, and the results 
             should therefore be interpreted with care."
           )
         ),
         modalDialog(title=topic,
                     "No help currently available for this topic.")
  )
}

.testDescription <- function(){
  tagList(
    tags$h3("Enrichment tests"),
    tags$p("Several target enrichment tests were benchmarked (results in the 
    benchmark tab), and are described below. The best overall tests were 
    selected as default for the app, and although other implemented tests are
    made available in the app's advanced options, their usage is not 
    recommended."),
    tags$p("The tests differ in the inputs they use, both in terms of the
    target annotation as well as of the type of signal in which enrichment is
    looked for. On the signal side, tests denoted as 'binary' compare features 
    (genes or transcripts) in a given set (e.g. your significantly downregulated
    genes) to those in a background set (i.e. over-representation analysis),
    whereas tests denoted as 'continuous' instead rely on a numeric input 
    signal, such as the mangitude or significance of changes in an input 
    differential expression analysis (by default, the tests use the sign of the
    foldchange multiplied by the -log10(FDR), which is well correlated to logFC
    for genes with low intra-group variability, and more robust than the latter).
    On the annotation side, tests can alsoeither use set membership (i.e. 
    whether or not a given feature is a predicted miRNA target) or numeric 
    values, such as the number of binding sites harbored by a given feature, 
    or a repression score (i.e. the extent to which a given feature is 
    predicted to be repressed by a miRNA)."),
    tags$ul(tags$li(tags$h4("Default tests"),
      tags$ul(tags$li(tags$b("siteoverlap")," (binary signal, set membership):",
                      tags$br(), "The siteoverlap test is based on Fisher's 
                      exact test, but using the number of sites on predicted 
                      targets and in the background instead of counting each
                      feature as one. While in theory this violates the 
                      assumption of independence of the counts (since all the
                      binding sites of a given transcript are either in or out
                      of the set), leading to slightly anti-conservative 
                      p-values, in practice this test is excellent at 
                      identifying the most enriched miRNA."),
              tags$li(tags$b("areamir")," (continuous signal, score or set 
                      membership):", tags$br(), 
                      "The areamir test is based on the analytic Rank-based 
                      Enrichment Analysis (aREA) test implemented in the 
                      'msviper' function of the ", tags$a("viper", target="_blank",
                      href="https://www.bioconductor.org/packages/devel/bioc/html/viper.html"),
                      " package. The test is akin to an analytical version of 
                      GSEA (see below), but it can additionally use degrees or 
                      likelihood of set membership. If repression scores are 
                      available in the annotation, areamir will therefore use a
                      (trimmed) version of it as set membership likelihood.")
      )),
      tags$li(tags$h5("Other tests implemented and evaluated"), tags$ul(
        tags$li(tags$b("overlap")," (binary signal, set membership):",
          tags$br(), "This test is based on Fisher's exact test, using the 
          number of features (i.e. transcripts/genes) among predicted targets 
          vs in the background (and therefore ignoring any site-based 
          information."),
        tags$li(tags$b("woverlap")," (binary signal, set membership):",
                tags$br(), "This test is like the above 'overlap' test, but
          corrects for UTR length using the Wallenius method, as implemented in
          the ", tags$a("goseq", target="_blank",
          href="https://www.bioconductor.org/packages/devel/bioc/html/goseq.html"),
          "package."),
        tags$li(tags$b("Mann-Whitney (KS)")," (continuous signal, set membership):",
                tags$br(), "This is the Mann-Whitney (also known as Wilcoxon)
                non-parametric test comparing targets and non-targets."),
        tags$li(tags$b("Kolmogorov-Smirnov (KS)")," (continuous signal, set membership):",
                tags$br(), "This is the Kolmogorov-Smirnov test comparing the 
                signal distribution of targets vs non-targets."),
        tags$li(tags$b("modscore")," (continuous signal, repression score):",
                tags$br(), "This is a linear regression testing the relationship
                between the input signal and the corresponding repression score
                predicted for a given miRNA."),
        tags$li(tags$b("modsites")," (continuous signal, number of sites):",
                tags$br(), "This is a linear regression testing the relationship
                between the input signal and the number of predicted binding 
                sites for a given miRNA, correcting for UTR length."),
        tags$li(tags$b("GSEA"), " (continuous signal, set membership)",
          tags$br(), "This test uses the multi-level fast GeneSet Enrichment 
          Analysis (GSEA) implemented in the ",tags$a("fgsea", target="_blank",
          href="https://www.bioconductor.org/packages/devel/bioc/html/fgsea.html"),
          "package, which is highly successful for Gene Ontology enrichment 
          analysis. In the context of our benchmark, however, it performed very 
          poorly."),
        tags$li(tags$b("regmir"), tags$br(), 
          "The regmir test uses lasso-regularized linear regression, which is 
          characterized by a high specificity and low sensitivity. The test will
          used binary or continuous inputs, as well as binary set membership or
          predicted repression score, depending on the availability of the input.
          Note that this test is very slow to run.")
        ))
      )
  )
}

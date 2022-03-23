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
      differential expression analysis (DEA) result you also have the 
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
      non-targets. In the presence of a real miRNA effect (that is gene expression
       changes were primarily caused by a microRNA), we expect to see a 
      gradual dose-response pattern (targets with stronger binding sites also 
      have stronger log-foldchanges), as is the case in this example.<br/><br/>
      Note that we could also have gotten to this CD plot by clicking on the
      corresponding point in the enrichment plot of the previous tab!",
      "This concludes our tutorial. While there are additionnal options here 
      and there, these are the app's main functionalities."
    )
  )
}

.getHelpModal <- function(topic){
  switch(topic,
         overview=modalDialog(title="enrichMiR workflow overview", easyClose=TRUE, size = "l",
                              tags$img(src="overview.png",height = 525, width = 750,
                                       style="display: block; margin-left: auto; margin-right: auto;"),
                              tags$head(tags$style("#shiny-modal img { max-width: 100%; }")),
                              tags$p("For further information regarding the individual steps consult the 
                                     specifi help pages or conduct the interactive tour."),
                              tags$p("This website is free and open to all users and there is no login requirement.")
         ),
         collections=modalDialog(title="Binding sites collections", easyClose=TRUE,
                                 tags$p(
"Each collection represents a set of binding sites (validated or predicted) for
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
               "miRNA binding site predictions using the ",
               tags$a(href="https://www.bioconductor.org/packages/devel/bioc/html/scanMiR.html",
                      "scanMiR", target="_blank"),
               "algorithm. ScanMiR employs a transcript-level affinity prediction method and is 
               compatible with both gene and transcript IDs."),
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
         deaformat=modalDialog(title="Required format of the DEA table", easyClose=TRUE,
           tags$p("The app should recognized the differential expression tables
                  generated from common RNAseq DEA package and saved as 
                  character-delimited files, assuming that the gene names/IDs 
                  are in the first column (which is normally the case).
                  In case of problem, try to format your table in this way:"),
           tags$pre(
"Gene,logFC,PValue,FDR
CERS2,-1.54,1.37e-19,8.85e-16
NR6A1,-1.76,3.92e-19,1.77e-15
TAGLN2,-1.92,4.47e-17,1.08e-13
FTL,-1.69,1.68e-15,2.05e-12
TPM4,-1.41,1.83e-10,6.86e-08
..."),
           tags$p("The character delimiter needs not be commas, and should be
           automatically recognized. With most target annotations, both Ensembl 
           IDs and gene symbols will be recognized, and the scanMiR annotation 
           additionally supports transcript-level DEAs with Ensembl transcript IDs."),
           tags$p("The character delimiter needs not be commas, and should be
           automatically recognized. With most target annotations, both Ensembl 
           IDs and gene symbols will be recognized, and the scanMiR annotation 
           additionally supports transcript-level DEAs with Ensembl transcript IDs."),
           tags$p("Because the accuracy of the fold-change estimates increases
             with the genes' read count, we recommend filtering the DEA table
             to include only highly expressed genes (e.g. the top 5000) before
             upload.")
         ),
         enrichplot=modalDialog(title="Enrichment plot", easyClose=TRUE, tags$p(
           "The enrichment plot shows the significance of the enrichment on the
           y-axis, and the magnitude of the enrichment on the x-axis. For basic
           overlap/over-representation tests, the enrichment is the fold 
           increase over the overlap expected by chance. For continuous tests
           (based on the results of a differential expression analysis), the 
           enrichment represents a magnitude of association with the 
           log-foldchanges of the putative targets."), 
           tags$p("In general, we recommend regarding the statistical significance 
                  of individual posttranscriptional regulators obtained via the 
                  different enrichment tests rather as rank instead of interpreting 
                  them as finite results. Hence high-ranking miRNAs/families are 
                  suggested to exert a functional role in the gene set at glance, 
                  though certainly need to be further validated."),
           tags$p("The size of the 
           points is proportional to the number of targets with binding sites for 
           the given miRNA/family. If miRNA expression data was given, the 
           points are colored by this expression."),
           tags$p("Hovering on the points/families will provide extra 
                  information about the family members and the test statistics."),
           tags$p("In case a DEA has been provided as input, double-clicking on a 
                  miRNA/family will automatically open the corresponding CD-plot."),
           tags$p("The plot size can be manually adapted by dragging the lower right part 
                  with the mouse.")
         ),
         browsercompatibility=modalDialog(title="Browser compatibility", easyClose=TRUE,
                        fluidRow(column(8, align="center", tableOutput("browsercomp")))
         ),
         cdplot=modalDialog(title="Cumulative distribution (CD) plots", easyClose=TRUE,
                            tags$p(
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
            happen when using TargetScan ", tags$em("conserved"), "sites (or
            miRTarBase sites), in which case you can try using all predicted 
            sites instead."),
           tags$p("Because the accuracy of the fold-change estimates increases
             with the genes' read count, we recommend using only highly 
             expressed genes to plot cummulative distributions of predicted 
             miRNA bindings, e.g. the top 5'000 - 7'000 expressed genes of a 
             dataset."),
           tags$p("Note that 'no site' indicates genes which are not annotated
             to have a binding site for that miRNA in the collection used. If
             the collection used is, for instance, the TargetScan conserved
             sites, there might nevertheless be non-conserved sites."),
           tags$p("The plot size can be manually adapted by dragging the lower 
              right part with the mouse.")
         ),
         tests=modalDialog(title="Enrichment tests", easyClose=TRUE, 
              .testIntro(), tags$p(
              "For more information on the tests, please see the documentation."
            )
         ),
        tests2=modalDialog(title="Enrichment tests", easyClose=TRUE, 
                           .testIntro(), .testDescription()
        ),
        testsadvanced=modalDialog(title="Enrichment tests (advanced)", easyClose=TRUE,
                                  .testDescription()),
        noExample=modalDialog(easyClose=TRUE,
          title="Example not available for this species", tags$p(
          "There is unfortunately no example data available for this species.
          You may however view an example by first switching to another species.")
        ),
        background=modalDialog(easyClose=TRUE, title="Choice of the background",
          tags$p("The background or universe defines the null hypothesis against
            is measured and tested an enrichment in the set of genes of 
            interest. As such, it applies only to tests that are based on 
            over-representation (e.g. siteoverlap), and not to tests that are 
            based on a continuous signal (e.g. areamir). In the context of
            over-representation analysis, the background is critical to 
            meaningful results. An inappropriate background, in particular a 
            background that is too broad (e.g. all genes), will often lead to 
            spurious results."),
          tags$p("The choice of the appropriate background depends on the 
            circumstances, but a key consideration is what genes ",tags$em(
            "could have made it"), " into your selection. For example, if your
            selection are genes that are differentially-expressed in a given 
            setting, then only genes that could have been differentially-
            expressed should be included in the background. This means that 
            genes which are too lowly-expressed for any differential pattern to
            be detected ought to be removed from the background.")
        ),
        modalDialog(title=topic, "No help currently available for this topic.")
  )
}

.testIntro <- function(){
  tagList(
    tags$p(
      "enrichMiR implements different statistical tests for target 
      enrichment. Some tests depend on continous inputs (e.g. fold-changes), 
      and hence require the results of a differential expression analysis (DEA)
      to be provided), while others are binary, i.e. based on the membership 
      in a gene set. When the input is a DEA, each binary tests will be 
      performed on both significantly upregulated and downregulated genes."
    ), tags$p(
      "Note that the tests were developed and benchmarked for
       application with predicted miRNA targets, not RNA binding proteins
       (RBPs), which have more degenerate binding patterns. Although we 
       offer the possibility to perform RBP target enrichment analyses, 
       the tests were not benchmarked for that purpose, and the results 
       should therefore be interpreted with care."
    ), tags$p(
      "Some tests gave poor performances in the benchmark, and should 
       therefore not be used. By default, the 'siteoverlap' and 'areamir' 
       tests are enabled -- these are the tests that gave the best performance.
      In addition, the binary version of the 'regmir' test provided excellent
      error control, albeit with a lower sensitivity."
    ) 
  )
}

.testDescription <- function(onlyDesc=FALSE){
  desc <- tagList(
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
    On the annotation side, tests can also either use set membership (i.e. 
    whether or not a given feature is a predicted miRNA target) or numeric 
    values, such as the number of binding sites harbored by a given feature, 
    or a repression score (i.e. the extent to which a given feature is 
    predicted to be repressed by a miRNA)."),
    tags$p("As is the case for over-representation analysis in general, for 
           tests based on binary inputs the choice of a good background is 
           critical. In many contexts, the background set of genes will be the 
           set of genes expressed in the system of interest."),
    tags$ul(tags$li(tags$h4("Default tests"),
                  tags$ul(
                    tags$li(tags$b("siteoverlap")," (binary signal, set membership):",
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
          information)."),
            tags$li(tags$b("woverlap")," (binary signal, set membership):",
                    tags$br(), "This test is like the above 'siteoverlap' test, 
          but corrects for UTR length using the Wallenius method, as implemented in
          the ", tags$a("goseq", target="_blank",
                        href="https://www.bioconductor.org/packages/devel/bioc/html/goseq.html"),
                    "package. The test performs similarly to the siteoverlap test."),
            tags$li(tags$b("Mann-Whitney (MW)")," (continuous signal, set membership):",
                    tags$br(), "This is the Mann-Whitney (also known as Wilcoxon)
                non-parametric test comparing targets and non-targets. This
                test performs badly in benchmarks and should not be used."),
            tags$li(tags$b("Kolmogorov-Smirnov (KS)")," (continuous signal, set membership):",
                    tags$br(), "This is the Kolmogorov-Smirnov test comparing the 
                signal distribution of targets vs non-targets.  This
                test performs badly in benchmarks and should not be used."),
            tags$li(tags$b("modscore")," (continuous signal, repression score):",
                    tags$br(), "This is a linear regression testing the relationship
                between the input signal and the corresponding repression score
                predicted for a given miRNA."),
            tags$li(tags$b("ebayes")," (continuous signal, repression score):",
                    tags$br(), "This is akin to the `modscore` tests, but performed
                using limma's moderated t-statistics."),
            tags$li(tags$b("lmadd")," (continuous signal, repression score):",
                  tags$br(), "This is the `ebayes` tests, followed by consecutive
                  fits adding each top miRNAs to the previous ones in a single 
                  model. This is especially useful to identify candidates which 
                  are not redundant with the top hit."),
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
                    "The regmir test uses constrained lasso-regularized regression,
          which has a high specificity but lower sensitivity. The test will use 
          binary or continuous inputs (using then either linear or binomial 
          regression), as well as binary set membership or predicted repression 
          score, depending on the availability of the input. The binary version
          of the test has shown the best performances.")
          ))
    )
  )
  if(onlyDesc) return(desc)
  tagList(
    tags$h3("Description of the enrichment tests"),
    tags$p("Several target enrichment tests were benchmarked (results in the 
    benchmark tab), and are described below. The best overall tests were 
    selected as default for the app, and although other implemented tests are
    made available in the app's advanced options, their usage is not 
    recommended."),
    desc
  )
}


.getTestsTable <- function(){
  data.frame(
    test=c("overlap","woverlap","siteoverlap","regmir (binary)","areamir","modSites",
           "ebayes","lmadd","regmir (continuous)","KS","MW","GSEA"),
    "input signal type"=rep(c("binary (set)", "continuous (DEA)"), c(4,8)),
    "annotation type"=c("genesets","nb sites","nb sites","genesets","scores",
                        "nb sites","scores","scores","scores",rep("genesets",3)),
    "description"=c(
      "Over-representation (ORA) of target genes among set",
      "ORA of binding sites, correcting for UTR length",
      "ORA of binding sites",
      "Regularized constrained logistic regression",
      "Score-weighted analytic rank enrichment analysis",
      "Linear regression of logFCs on nb of binding sites",
      "Linear regression of logFCs on predicted repression scores",
      "Consecutive additive linear regression of logFCs on predicted repression scores",
      "Regularized constrained lienar regression of logFCs on predicted repression scores",
      "Kolmogorov-Smirnov (KS) test on logFCs",
      "Mann-Whitney / Wilcoxon test on logFCs",
      "Gene set enrichment analysis (GSEA)"))
}


.getBrowserCompTable <- function(){
  data.frame(
    "OS"=c("Linux","MacOS","Windows"), 
    "Version" = c("Ubuntu","Big Sur","10"),
    "Chrome" = c("96.0","96.0","89.0"),
    "Firefox" = c("95.0","95.0","86.0"),
    "Microsoft Edge" = c("n/a","n/a","89.0"),
    "Safari" = c("n/a","15.0","n/a")
     )
}


.exampleBackground <- function(){
  c("PLEKHB2", "FBXO21", "PIGS", "SLC7A1", "YKT6", "RABL6", "SLC1A5", 
    "OCRL", "SH3RF1", "SLC9A1", "SLC7A5", "SRPRA", "G6PC3", "ARL2", 
    "CTDSPL", "CSRP1", "TBC1D22B", "ZER1", "SLC52A2", "GNPDA1", "PRDM4", 
    "VAMP3", "TMEM87A", "SEPTIN2", "SHCBP1", "RELL1", "PRTFDC1", 
    "CPEB1", "PHKA1", "RNF38", "TMED8", "CASP7", "PROSER1", "LMNB2", 
    "DSTYK", "TMEM216", "RAVER1", "PRUNE1", "ASF1B", "NCDN", "TGFBRAP1", 
    "CS", "MTHFD2", "SFT2D1", "ARHGAP1", "IQGAP1", "ATN1", "CTDNEP1", 
    "VPS37B", "PLCD3", "PKM", "POLR3D", "SLC25A6", "PRKRA", "HECTD3", 
    "SULF1", "SERINC5", "DYNC1LI2", "BCAT2", "VPS4A", "AL356776.1", 
    "DMAC1", "MROH8", "MSH2", "PANX2", "PAQR9", "SWI5", "DHX9P1", 
    "RAP1GDS1", "OLFM2", "PGF", "TBC1D12", "TNFSF12", "ANP32BP1", 
    "OSGIN2", "LMCD1", "CNOT9", "TTC6", "YTHDF2", "BOD1", "ZMYM3", 
    "USF1", "STK4", "PYCR3", "AC104083.1", "NSUN5P2", "PTAFR", "CFAP58", 
    "CDKN2AIP", "ARPC2", "HMGN2P46", "PJA1", "HSPA6", "AC079416.3", 
    "MYH9", "MRTFB", "SDK1", "HNRNPA1P69", "RAB10", "BUD13", "NDUFA2", 
    "STK31", "SEC11B", "ROMO1", "KDELC1P1", "EXOSC5", "AC092807.3", 
    "AC004492.1", "PRSS53", "MTRNR2L6", "RPS2P6", "RPL5P1", "AC013391.2", 
    "PNPT1P1", "CALM2", "AC027801.2", "ANO8", "RAB2B", "AL136038.4", 
    "PPM1B", "BHLHE41", "SBF1", "LIFR-AS1", "PSMB4", "POGK", "CSDC2", 
    "PPIAP26", "PHRF1", "NPM1P9", "PDIA2", "BCRP2", "GAPDHP58", "CLEC16A", 
    "AC234775.3", "EP400", "VTI1B", "PTPN6", "SIMC1", "AL034397.1", 
    "RBBP8", "PEX14", "ABCA11P", "TEN1-CDK3", "ZYG11B", "AGO1", "TAF3", 
    "GBX2", "PLIN5", "PKDCC", "SMNDC1", "RUSC2", "FAM20C", "HBP1", 
    "SDHAF1", "ZBTB40-IT1", "TMEM14C", "TUFM", "SH3GL2", "MDM4", 
    "MIR3936HG", "FAM210A", "SNRPGP10", "CCAR1", "RORB", "AL645940.1", 
    "RPN1", "LSM6", "PPRC1", "CD70", "AC027307.3", "ZNF789", "BRSK2", 
    "MEGF6", "DPYSL3", "CRYBG3", "DMAP1", "TMEM101", "CDKN1B", "UCP1", 
    "ZNF222", "LINC00926", "RGS7", "EML4", "PCBP1", "BTF3P8", "INTS1", 
    "AC106786.1", "RNF168", "CASP9", "MLXIPL", "CDX2", "ANKRD54", 
    "AC027682.7", "PRRG1", "AC209007.1", "AC026356.1", "AC104964.3", 
    "THG1L", "AC139530.1", "TEX22", "RPL7P59", "CRB2", "PDCD2", "HTR1D", 
    "PPIA", "YARS2", "AC092183.1", "TVP23C", "EEF1A1P6", "DNAAF1", 
    "INAFM1", "TRIOBP", "PIR", "NOP10", "DDX50P2", "ACTL8", "UQCRC2", 
    "RAB7A", "FBXL19-AS1", "PA2G4P4", "MAP4K2", "HADHAP2", "CKAP5", 
    "ZNF215", "GPX4", "MTX1", "TTLL4", "AHNAK", "KRT18P31", "ELP4", 
    "NCAM1", "DOK4", "ACAA2", "AC092431.2", "RPL22P8", "SP9", "ST6GALNAC6", 
    "TRA2A", "NOP2", "PLXNB2", "GADD45GIP1", "ASH1L-AS1", "ZNF843", 
    "RPS23P8", "KARS", "ZNF407", "TRAPPC2B", "ZNF518A", "KCNA3", 
    "UCKL1", "NPM1P6", "NOL3", "SPR", "NOP14", "LGMN", "VWA5A", "GAPDHP25", 
    "ARTN", "KBTBD7", "DCXR", "ANO6", "AC008147.4", "HAX1", "NEDD1", 
    "ZNF772", "AC020910.5", "DEPDC7", "HERC2P2", "DCHS1", "RPL37P1", 
    "PLIN1", "APLF", "ASB13", "NOP56", "PCSK4", "AMACR", "PLCL1", 
    "G0S2", "VASN", "PPP3CB", "PSMD4", "AL512791.2", "STK17A", "TSSK3", 
    "AC006460.2", "PPP1R26", "ELF1", "ACTBP8", "CSNK2B", "EDRF1", 
    "AC027097.1", "PTBP1", "KLLN", "GSTK1", "VKORC1L1", "COL11A2", 
    "C15orf62", "PRKCG", "ACVR2B-AS1", "RDH16", "WASHC2C", "AC114980.1", 
    "TXN2", "RPL19P16", "SRR", "RPS19BP1", "AC010680.5", "RN7SL832P", 
    "PCNA", "AC013403.2", "UBE2Q2", "AC016737.1", "TMPO", "API5P2", 
    "FAAH2", "STAMBP", "AL445487.1", "RPS2P44", "SRRM2-AS1", "AC079193.2", 
    "GATAD2A", "FEZF1", "HHEX", "RPF1", "IRX3", "AC093788.1", "STON2", 
    "TCN2", "CMSS1", "RPL27AP5", "DPH1", "AL136116.3", "TOMM5", "RPL14", 
    "AL365295.1", "AC090015.1", "VPS13B", "MXD4", "CHP2", "RN7SL535P", 
    "ACTG1P14", "AC087343.1", "PRPF38AP2", "GRIPAP1", "AC107068.1", 
    "CTSL", "CTSF", "SNHG11", "HHIP-AS1", "NMUR1", "R3HCC1L", "PPP1R13L", 
    "UGDH-AS1", "PDF", "GINS4", "PIPSL", "AC109347.1", "SLC16A10", 
    "RPS18P13", "ASPHD2", "EEF1E1P1", "AL096803.2", "RIC8B", "TFDP1", 
    "RBM6", "METTL4", "SULT4A1", "RPS2P46", "DSTN", "AC005072.1", 
    "ACAA1", "AC012313.1", "GSDME", "HNRNPA1P54", "TOGARAM2", "AC022498.2", 
    "TVP23B", "RPS12P26", "AC131235.1", "RCCD1", "INPP1", "WDFY1", 
    "AC026367.2", "MCF2L-AS1", "SETP21", "TMEM138", "SNORA33", "PABPN1", 
    "OGDH", "AL122023.1", "CHEK1", "CCN2", "AL139260.1", "SYMPK", 
    "CDKN2AIPNL", "MBD5", "NCOA4P4", "BLVRB", "DSE", "MIEF2", "RBM8B", 
    "ZNF350", "ADSL", "PARN", "LNCTAM34A", "PNMA6A", "PIGW", "TIMM29", 
    "IL21R-AS1", "NBPF3", "CYB5RL", "CERT1", "MAP7D2", "RUNX1", "AC127024.8", 
    "JAKMIP2", "NIBAN2", "TIGAR", "RNF227", "UNG", "MTOR", "TMEM107", 
    "JUND", "GEMIN7", "CALCOCO2", "TOR2A", "ATAD3A", "C1orf43", "CRK", 
    "AL353151.1", "GAPDHP23", "DEGS1", "KXD1", "AF131216.4", "GTF2H2C", 
    "LAMP2", "IPO7P1", "RAP2C-AS1", "LFNG", "MTND2P9", "CNOT7", "AC134349.1", 
    "GMPSP1", "NAA25", "MMP24OS", "VMP1", "NTNG2", "MCCC1", "RPS25", 
    "ILF2P2", "AL035530.2", "AC024451.1", "GNMT", "XRCC6P1", "NSMCE3", 
    "AC026124.1", "AC074194.1", "NHP2P1", "CHD2", "AC073046.1")
}

.exampleGeneset <- function(){
  head(.exampleBackground(),60)
}

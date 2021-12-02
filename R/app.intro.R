.getAppIntro <- function(){
  data.frame(
    element=c("#menu_species", "#collection_input", "#exprMirs_box", 
              "#expressed_miRNAs_box", "#menu_input", "#input_type",
              "#menu_input", "#menu_enrich", "#enrich"),
    intro=c(
      "The first step in using the enrichMiR app is to make sure that you're
       working with the right species and annotation. To set this, you open 
       the 'Species and miRNAs' tab.",
      "Once you've selected a species, you will be able to see and select one
       of its available target annotations (see the help button for information
       about the annotations). For the sake of this example, select the human
       'targetScan conserved miRNA BS'.",
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
      #### add btn to get example DEA
      "If you've provided a valid input, the input menu on the left will show
      a blue badge instead of a red one.",
      "We're now ready to launch the enrichment analysis!<br/>First go to the 
      'Enrichment analysis' tab",
      "Then click the 'Run enrichMir!' button. The analysis will take a moment."
    )
  )
}
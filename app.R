# By Avery Davis Bell
# This is a Shiny web application. You can run the application locally by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Purpose: plot gene counts from DESeq2 object of CeNDR data

require(shiny, quietly = T)
require(ggplot2, quietly = T)
require(data.table)
require(DT)
require(formattable)

#### Functions ####
# --- Data, input, etc
getgenename<-function(ingene, geneinfo){
  # Determines gene_id, display_name, sequence_name, and what to show on app from ingene, which might correspond to any of those gene identifiers
  # In: ingene, gene_id (WBGene...), locus name (e.g. pos-1), or sequence_name (like ZK...) for which to get full name info
  #     geneinfo, data.table with one row per gene. Columns gene_id, display_name (pos-1 like name unless gene doesn't have one, in which case it's sequence name), sequence_name
  # Out: one-row data.table with gene information. Columns:
  #   gene_id, display_name, sequence_name - as in geneinfo input
  #   app_info - gene name/info to display to user: display name and parenthetically gene ID; optionally also sequence name if sequence name is provided and different from display name
  
  if(geneinfo[,ingene%in%gene_id]){
    out<-geneinfo[gene_id==ingene, .(gene_id, display_name, sequence_name,
                                     app_info = paste0(display_name, " (", gene_id, ")"))]
  }else if(geneinfo[,ingene%in%display_name]){
    out<-geneinfo[display_name==ingene, ][1, ] # currently display_names are unique but there's a possibility they wouldn't be
    out[,app_info:=paste0(display_name, " (", gene_id, ")")]
  }else if(geneinfo[,ingene%in%sequence_name]){ # Not in display_name, so want to show both in app_info
    out<-geneinfo[sequence_name==ingene, ][1, ] # a few sequence_names map to multiple genes - though not currently if they don't have display_name, still worth catching
    out[,app_info:=paste0(display_name, "/", sequence_name, " (", gene_id, ")")]
  }else{
    
    out<-NULL # consider change to whole DT with NAs?. Gene doesn't have expression data/ not analyzed here [wrong classification etc?].
  }
  return(out)
}

# --- Plotting
txtplot<-function(mytxt, mycex = 1.5, mycolor = "red", mytitle = ""){
  # Makes blank plot with text mytxt displayed in middle. One way to consider conveying info.
  # In: mytxt, text to plot
  #     mycex, cex (size factor) for text
  #     mycolor, color for text
  #     mytitle, optional title for plot
  plot(1:10, 1:10, type = 'n', axes = F, xlab = "", ylab = "", main = mytitle)
  text(5, 5, mytxt, cex = mycex, col = mycolor)
}

exprhist<-function(exprinfo, mygene, genetitle = "", mysubt = ""){
  # Plots MEAN expression of hist across strains
  # In: exprinfo, full data
  #     mygene, gene_id for gene to plot
  #     genetitle, title for plot
  #     mysubt, subtitle for plot
  pdat<-unique(exprinfo[gene_id==mygene, .(strain, meanexp)])
  
  plt<-ggplot(pdat, aes(meanexp)) + geom_histogram() +
    xlab("Expression (mean in strain; normalized; log2)") + 
    ylab("Number of strains") +
    ggtitle(genetitle, subtitle = mysubt) + theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 15), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          legend.title = element_text(size=15), legend.text = element_text(size=14),
          strip.text.x = element_text(size = 15), title = element_text(size = 17),
          legend.position = "none")
    
  
  return(plt)
}

allstrncounts<-function(exprinfo, onestrain, mygene, genetitle = "", mysubt = ""){
  # Plots focal strain & strains with highest and lowest mean expression as dot plot!!
  
  # Data fussing
  gdat<-unique(exprinfo[gene_id==mygene, .(strain, meanexp)])
  plotstrains<-c(gdat[gdat[,which.min(meanexp)], strain], onestrain, gdat[gdat[,which.max(meanexp)], strain])
  pdat<-exprinfo[gene_id==mygene & strain%in%plotstrains, .(strain, sampexp)]
  
  # Plot
  plt<-ggplot(pdat, aes(strain, sampexp)) + geom_dotplot(aes(fill = strain), binaxis = "y", stackdir = "center") +
    xlab("Strain") + ylab("log2(normalized expression)") + ggtitle(genetitle, subtitle = mysubt) + 
    theme_bw() +
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 15), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 14),
          legend.title = element_text(size=15), legend.text = element_text(size=14),
          strip.text.x = element_text(size = 15), title = element_text(size = 17),
          legend.position = "none")
  
  # ?? Un-log numbers on Y axis (keep same scale, just label with actual values?)
  return(plt)
}

miscgenestraininfo<-function(exprinfo, onestrain, mygene){
  # Generates text describing various aspects of this gene
  
  # Get values
  thisdat<-unique(exprinfo[gene_id==mygene, 
                           .(strain, gene_id, display_name, biotype, chr, start, end, strand, hypdiv.strain, hypdiv.any, meanexp)])
  mymean<-round(thisdat[strain==onestrain, meanexp], 2)
  myrank<-thisdat[order(meanexp), which(strain==onestrain)]
  mingene<-round(thisdat[,min(meanexp)], 2)
  maxgene<-round(thisdat[,max(meanexp)], 2)
  
  # Generate text
  outtxt<-paste0("<p>",
                c(paste0("This gene is located at ", thisdat[1, chr], ":", thisdat[1, start],"-", thisdat[1, end], "."),
                paste0("Is this gene hypervariable (previously termed hyperdivergent) in ", onestrain, "? ", thisdat[strain==onestrain, hypdiv.strain], ". ",
                       "Is this gene hypervariable (previously termed hyperdivergent) in any CeNDR 20220216 isotype? ", thisdat[strain==onestrain, hypdiv.any], "."),
                paste0("At this gene, minimum expression in any strain is ", mingene, ", and maximum expression in any strain is ", maxgene, "."),
                paste0("In strain ", onestrain, ", mean expression is ", mymean, ". This strain ranks ", myrank, "th of ", nrow(thisdat), " strains in expression.")),
                "</p>")
  
  # Return text
  return(outtxt)
}

#### Load data ####
load("cendr_rna_app_data.RData") # One big data.table, not too quick to load; titled exprinfo
strains.exp<-exprinfo[,unique(strain)]

#### UI ####
ui<-fluidPage(
  titlePanel("Gene expression for 208 wild C. elegans strains from CaeNDR", 
             windowTitle = "CaeNDR expression data"),
  fluidRow(column(12, h5("RNA-seq data from: Zhang et al 2022[1]"), 
                  h5("RNA quantifications (via strain-specific salmon and DESeq2 tximport/vst normalization) from: Bell and Paaby 2024[2]"))),
  
  sidebarLayout(
    sidebarPanel(
      h4("Select data to display"),
      textInput(inputId="usergene", label="Gene (WBGene ID, common name, or sequence name):", value="vit-6"),
      selectInput(inputId = "strn", label = "Specific strain of interest:", 
                  choices = strains.exp, selected = "CB4856"),
      # --- displayu selected gene info
      h5(textOutput({"geneformat"})),
      
      # --- Style of sidebar
      style = "position:fixed;width:30%;"
    ),
    
    mainPanel(
      # --- Expression histogram
      fluidRow(column(12), h3("Expression across all RNA-sequenced strains")),
      fluidRow(column(12), plotOutput(outputId = "exphistplt")),
      hr(),
      
      # --- count dot plot & explanation
      fluidRow(column(12), h3("Per-sample expression in selected strains")),
      fluidRow(column(12), "Strains shown are your selected strain, the strain with the lowest mean expression at this gene, and the strain with the highest mean expression at this gene."),
      fluidRow(column(12), plotOutput(outputId = "countdotplot")),
      hr(),
      
      # --- Extra gene info (text)
      fluidRow(column(12), h3("Other information about this gene (expression and non-expression based)")),
      fluidRow(column(12), htmlOutput(outputId = "otherinfo")),
      hr(),
      
      # --- citations
      fluidRow(column(12, h3("Citations and data sources"),
                      DT::dataTableOutput(outputId = "citations"))),
      hr(),
      
      # --- license etc
      fluidRow(column(12, h4("About"),
                      "This (",
                      tags$a(href = "https://shiny.rstudio.com/", "Shiny", target = "_blank", .noWS = "outside"), ") app written - and hopefully maintained - by ", 
                      tags$a(href = "https://www.averydavisbell.com", "Avery Davis Bell", target = "_blank", .noWS = "outside"), 
                      ", postdoc in ", 
                      tags$a(href = "https://genaamics.org", "Annalise Paaby's lab", target = "_blank", .noWS = "outside"), ". 2024."),
               column(12, tags$i("This work is licensed under a ",
                                 tags$a(href = "https://creativecommons.org/licenses/by/4.0/", "Creative Commons Attribution 4.0 International License", target = "_blank", .noWS = "outside"), "."))),
    )
  )
)

#### Server ####
server<-function(input, output){
  # --- Get gene name nicely formatted for display and use
  geneInput<-reactive(getgenename(ingene = input$usergene, geneinfo = unique(exprinfo[strain==input$strn, .(gene_id, display_name, sequence_name)])))
  
  output$geneformat<-renderText({
    if(is.null(geneInput())){ # gene not in input - display so
      "******The provided gene is not included in this dataset******"
    }else{
      paste("Selected gene:", geneInput()$app_info)
    }
  })
  
  output$exphistplt<-renderPlot({
    if(!is.null(geneInput())){
      exprhist(exprinfo, mygene = geneInput()$gene_id, genetitle = geneInput()$app_info, mysubt = "")
    }else{
      txtplot("Gene not included in dataset")
    }
  })
  
  output$countdotplot<-renderPlot({
    if(!is.null(geneInput())){
      allstrncounts(exprinfo, onestrain = input$strn, mygene = geneInput()$gene_id, genetitle = geneInput()$app_info, mysubt = "")
    }else{
      txtplot("Gene not included in dataset")
    }
  })
  
  output$otherinfo<-renderText({
    if(!is.null(geneInput())){
      miscgenestraininfo(exprinfo, onestrain = input$strn, mygene = geneInput()$gene_id)
    }else{
      "Gene not in dataset"
    }
  })
  
  # ---- Citations etc
  output$citations<-DT::renderDataTable({
    cits<-data.table(ref = c(1, 2),
                     Citation = c("<p>Zhang, G., N. M. Roberto, D. Lee, S. R. Hahnel and E. C. Andersen, 2022 The impact of species-wide gene expression variation on Caenorhabditis elegans complex traits. Nat Commun 13: 3462.</p>", 
                                  "<p>soon to be available at microPublication Biology</p>"),
                     Link = c("https://doi.org/10.1038/s41467-022-31208-4", "https://wildworm.biosci.gatech.edu/"),# link back to data page temporarily
                     seealso = c("https://caendr.org/", "https://wildworm.biosci.gatech.edu/")) 
    ## format links
    cits[, Link:=paste0("<a href='", Link, "' target='_blank'>", Link, "</a>")]
    cits[, Citation:=paste(Citation, Link)]
    cits[, Link:=NULL]
    cits[, seealso:=paste0("<a href='", seealso, "' target='_blank'>", c("CaeNDR", "Paaby lab interactive data"), "</a>")]
    ## format names
    setnames(cits, c("ref", "seealso"), c("Reference number", "See also"))
    as.datatable(formattable(cits), options = list(pageLength = nrow(cits), searching = F, paging = F, info = F, ordering = F),
                 rownames = F)
  })
}

#### Run the application ####
shinyApp(ui = ui, server = server)
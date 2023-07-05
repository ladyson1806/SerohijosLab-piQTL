library(shiny)
library(R.utils)
library(glue)
library(igvShiny)
library(BiocManager)
options(repos = BiocManager::repositories())

if(!dir.exists("tracks"))
  dir.create("tracks")
addResourcePath("tracks", "tracks")

PPI_table = read.csv('./data/PPI_reference_barcodes.csv')

# Define UI for app t ----
ui <- fluidPage(
  # Sidebar layout with input and output definitions ----
  # sidebarLayout(
  # Sidebar panel for inputs ----
  sidebarPanel(height=6,
               # Input: Slider for the number of bins ----
               column(6,
               selectInput(inputId = "ppi", h4("PPI of interest"),
                           choices = PPI_table$PPI)
               ),
               column(6,
                      selectInput(inputId = "drug", h4("Drug of interest"),
                                  choices = c("Fluconazole", "5.FC"," Metformin", "Trifluoperazine"))
               ),
               width=12,
               h4("Loading annotations"),
               actionButton("addTracks_0", "Add LD blocks, CUT, SUT and XUT annotations"),
               br(), br(),
               actionButton("addGwasTrackButton_DrugMTX", "Add piQTL Track (DRUG & MTX+)"),
               actionButton("addGwasTrackButton_noDrugMTX", "Add piQTL Track (noDRUG & MTX+)"),
               br(), br(),
               h4("Region of interest"),
               textInput("roi", label="", placeholder="Gene or chrN:start-end"),
               actionButton("searchButton", "Search")
               
               ),
    
    # Main panel for displaying outputs ----
    mainPanel(width=12,
              igvShinyOutput('igvShiny_0')
  )
)

# Server side ----
server <- function(input, output, session) {

    observeEvent(input$searchButton, {
    searchString = isolate(input$roi)
    printf("--- search: %s", searchString)
    if(nchar(searchString) > 0)
      showGenomicRegion(session, id="igvShiny_0", searchString)
  })
  
  output$igvShiny_0 <- renderIgvShiny({
    cat("--- starting renderIgvShiny\n");
    printf("---- addGFF3Track")
    genomeOptions <- parseAndValidateGenomeSpec(genomeName="sacCer3")
    x <- igvShiny(genomeOptions,
                  tracks=list(LDs)
    ) 

    f_LDs <-"./data/genome_annotations/LD_blocks.csv"
    LDs.gff3 <- read.table(f_LDs, sep=",", header=TRUE, as.is=TRUE)
    LDTrack <- loadGFF3TrackFromLocalData(session, id="igvShiny_0", "LD blocks", LDs.gff3, color='grey', colorByAttribute='type', colorTable=list(), displayMode="EXPANDED", trackHeight=50, visibilityWindow=100000)
    display(LDTrack, session, id="igvShiny_0", deleteTracksOfSameName = TRUE)
    
    f_CUTs <- "./data/genome_annotations/Xu_2009_CUTs_V64.csv"
    CUTs.gff3 <- read.table(f_CUTs, sep=",", header=TRUE, as.is=TRUE)
    loadGFF3TrackFromLocalData(session, id="igvShiny_0", "CUTs Xu 2009", CUTs.gff3, color='darkgreen', colorByAttribute='type', colorTable=list(), displayMode="EXPANDED", trackHeight=50, visibilityWindow=100000)
    
    f_SUTs <- "./data/genome_annotations/Xu_2009_SUTs_V64.csv"
    SUTs.gff3 <- read.table(f_SUTs, sep=",", header=TRUE, as.is=TRUE)
    loadGFF3TrackFromLocalData(session, id="igvShiny_0", "SUTs Xu 2009", SUTs.gff3, color='green', colorByAttribute='type', colorTable=list(), displayMode="EXPANDED", trackHeight=50, visibilityWindow=100000)
    
    f_XUTs <- "./data/genome_annotations/van_Dijk_2011_XUTs_V64.csv"
    XUTs.gff3 <- read.table(f_XUTs, sep=",", header=TRUE, as.is=TRUE)
    loadGFF3TrackFromLocalData(session, id="igvShiny_0", "XUTs Van Dijk 2011", XUTs.gff3, color='palegreen', colorByAttribute='type', colorTable=list(), displayMode="EXPANDED", trackHeight=50, visibilityWindow=100000)
    
    cat("--- ending renderIgvShiny\n");
    return(x)
  })
  
  observeEvent(input$addTracks_0, {
    # indexURL <- "https://s3.amazonaws.com/igv.org.genomes/SacCer3/Saccharomyces_cerevisiae.R64.109.gff3.gz.tbi"
    printf("---- addGFF3Track")
    f_LDs <-"./data/genome_annotations/LD_blocks.csv"
    LDs.gff3 <- read.table(f_LDs, sep=",", header=TRUE, as.is=TRUE)
    loadGFF3TrackFromLocalData(session, id="igvShiny_0", "LD blocks", LDs.gff3, color='grey', colorByAttribute='type', colorTable=list(), displayMode="EXPANDED", trackHeight=50, visibilityWindow=100000)
    
    
    f_CUTs <- "./data/genome_annotations/Xu_2009_CUTs_V64.csv"
    CUTs.gff3 <- read.table(f_CUTs, sep=",", header=TRUE, as.is=TRUE)
    loadGFF3TrackFromLocalData(session, id="igvShiny_0", "CUTs Xu 2009", CUTs.gff3, color='darkgreen', colorByAttribute='type', colorTable=list(), displayMode="EXPANDED", trackHeight=50, visibilityWindow=100000)
    
    f_SUTs <- "./data/genome_annotations/Xu_2009_SUTs_V64.csv"
    SUTs.gff3 <- read.table(f_SUTs, sep=",", header=TRUE, as.is=TRUE)
    loadGFF3TrackFromLocalData(session, id="igvShiny_0", "SUTs Xu 2009", SUTs.gff3, color='green', colorByAttribute='type', colorTable=list(), displayMode="EXPANDED", trackHeight=50, visibilityWindow=100000)
    
    f_XUTs <- "./data/genome_annotations/van_Dijk_2011_XUTs_V64.csv"
    XUTs.gff3 <- read.table(f_XUTs, sep=",", header=TRUE, as.is=TRUE)
    loadGFF3TrackFromLocalData(session, id="igvShiny_0", "XUTs Van Dijk 2011", XUTs.gff3, color='palegreen', colorByAttribute='type', colorTable=list(), displayMode="EXPANDED", trackHeight=50, visibilityWindow=100000)
    
  })
    
  observeEvent(input$addGwasTrackButton_DrugMTX, {
    printf(glue("---- Adding GWASTrack for {input$ppi} with {input$drug} (MTX+)"))
    # colors <- as.list(rep(c("black", "lightgray"), 17))
    f <- glue("https://raw.githubusercontent.com/ladyson1806/public_hosting/main/piQTL_mapping/formatted_for_genome_browser/{input$ppi}_MTX_{input$drug}_avg_logratio_Fitness_minus_ref.csv")
    tbl.gwas <- read.csv(f, sep='\t')
    tbl.gwasTrack <- GWASTrack(glue('{input$ppi} with {input$drug} (MTX+)'), tbl.gwas, chrom.col=2, pos.col=3, pval.col=6, trackHeight=200)
    display(tbl.gwasTrack, session, id="igvShiny_0", deleteTracksOfSameName = TRUE)
  })
  
  observeEvent(input$addGwasTrackButton_noDrugMTX, {
    printf(glue("---- Adding GWASTrack for {input$ppi} without {input$drug} (MTX+)"))
    # colors <- as.list(rep(c("black", "lightgray"), 17))
    f <- glue("https://raw.githubusercontent.com/ladyson1806/public_hosting/main/piQTL_mapping/formatted_for_genome_browser/{input$ppi}_MTX_noDrug_avg_logratio_Fitness_minus_ref.csv")
    tbl.gwas <- read.csv(f, sep='\t')
    tbl.gwasTrack <- GWASTrack(glue('{input$ppi} without {input$drug} (MTX+)'), tbl.gwas, chrom.col=2, pos.col=3, pval.col=6, trackHeight=200)
    display(tbl.gwasTrack, session, id="igvShiny_0", deleteTracksOfSameName = TRUE)
  })

}
# Create Shiny app ----

  
shinyApp(ui = ui, server = server)

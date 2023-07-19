library(shiny)
library(glue)
library(BiocManager)
library(igvR)

options(repos = BiocManager::repositories())
options(shiny.host = '0.0.0.0')
# options(shiny.port = 3839)


# Define UI for app that draws a histogram ----

PPI_table = read.csv('./data/PPI_reference_barcodes.csv')

ui <- fluidPage(
  
  # App title ----
  titlePanel("piQTL - Interactive Genome Browser"),
  
  # Sidebar layout with input and output definitions ----
  # sidebarLayout(
  # Sidebar panel for inputs ----
  sidebarPanel(width=12,
               # Input: Slider for the number of bins ----
               selectInput(inputId = "ppi", h4("PPI of interest"),
                           choices = PPI_table$PPI),
               selectInput(inputId = "drug", h4("Drug of interest"),
                           choices = c("5.FC", "Fluconazole", "Metformin", "Trifluoperazine"))
               ),
    
    # Main panel for displaying outputs ----
    mainPanel(width=12,
    plotOutput("igv_viewer")
  )
)


romanize <- function(x) {
  if (class(x) == 'integer') as.roman(x)
  else x
}

server <- function(input, output) {
  MTX_GWAS_table <- reactive({
    PPI <- input$ppi
    DRUG <- input$drug
    MTX_infile = glue('./data/piQTLs/piQTL_without_MTX/{PPI}_MTX_{DRUG}_avg_logratio_Fitness_minus_ref.csv')
    MTX_table <- read.csv(MTX_infile)
    MTX_table['CHR_igv'] <- lapply(MTX_table['CHR'], romanize)
    MTX_table = MTX_table[c("SNP", "CHR_igv", "BP", "EFFECTSIZE", "SE", "P", "snps_class_up", "snps_class_down", "locus_id", "GENE", "sgd_id")]
    MTX_table
  })
  
  
  noDrug_MTX_GWAS_table <- reactive({
    PPI <- input$ppi
    DRUG <- input$drug
    noDrug_MTX_infile = glue('./data/piQTLs/piQTL_without_MTX/{PPI}_MTX_noDrug_avg_logratio_Fitness_minus_ref.csv')
    noDrug_MTX_table <- read.csv(noDrug_MTX_infile)
    noDrug_MTX_table['CHR_igv'] <- lapply(noDrug_MTX_table['CHR'], romanize)
    noDrug_MTX_table = noDrug_MTX_table[c("SNP", "CHR_igv", "BP", "EFFECTSIZE", "SE", "P", "snps_class_up", "snps_class_down", "locus_id", "GENE", "sgd_id")]
    noDrug_MTX_table
  })
  
  noDrug_noMTX_GWAS_table <- reactive({
    PPI <- input$ppi
    DRUG <- input$drug
    noDrug_noMTX_infile = glue('./data/piQTLs/piQTL_without_MTX/{PPI}_noMTX_noDrug_avg_logratio_Fitness_minus_ref.csv')
    noDrug_noMTX_table <- read.csv(noDrug_noMTX_infile)
    noDrug_noMTX_table['CHR_igv'] <- lapply(noDrug_noMTX_table['CHR'], romanize)
    noDrug_noMTX_table = noDrug_noMTX_table[c("SNP", "CHR_igv", "BP", "EFFECTSIZE", "SE", "P", "snps_class_up", "snps_class_down", "locus_id", "GENE", "sgd_id")]
    noDrug_noMTX_table
  })

  # noMTX_GWAS_table <- reactive({
  #   PPI <- input$ppi
  #   DRUG <- input$drug
  #   noMTX_infile = glue('../../data/QTL/{PPI}_noMTX_{DRUG}_avg_logratio_Fitness_minus_ref.csv')
  #   noMTX_table <- read.csv(noMTX_infile)
  #   noMTX_table['CHR_igv'] <- lapply(noMTX_table['CHR'], romanize)
  #   noMTX_table = noMTX_table[c("SNP", "CHR_igv", "BP", "EFFECTSIZE", "SE", "P", "snps_class_up", "snps_class_down", "locus_id", "GENE", "sgd_id")]
  #   noMTX_table
  # })
  
  igv <- igvR(host = '0.0.0.0')
  setBrowserWindowTitle(igv, "piQTL - Interactive Genome Browser")
  setGenome(igv, "sacCer3")
  
  output$igv_viewer <- renderPlot({
    
    
    LD.gff3 <- read.table("./data/genome_annotations/LD_blocks.csv",
                            sep=",", as.is=TRUE)
    colnames(LD.gff3) <- c("seqid", "source", "type", "start", "end", "score", "strand",
                            "phase", "attributes")
    LD_track <- GFF3Track("LD blocks 2023", LD.gff3,
                           url=NA_character_, indexURL=NA_character_, displayMode="EXPANDED", trackHeight=70,
                           visibilityWindow=100000)
    displayTrack(igv, LD_track)
    
    
    XUTs.gff3 <- read.table("./data/genome_annotations/van_Dijk_2011_XUTs_V64.csv",
                           sep=",", as.is=TRUE)
    XUT_track <- GFF3Track("XUTs Van Dijk 2011", XUTs.gff3,
                       url=NA_character_, indexURL=NA_character_, displayMode="EXPANDED", trackHeight=70,
                       visibilityWindow=100000)
    displayTrack(igv, XUT_track)
    
    
    CUTs.gff3 <- read.table("./data/genome_annotations/Xu_2009_CUTs_V64.csv",
                            sep=",", as.is=TRUE)
    CUT_track <- GFF3Track("CUTs Xu 2009", CUTs.gff3,
                           url=NA_character_, indexURL=NA_character_, displayMode="EXPANDED", trackHeight=70,
                           visibilityWindow=100000)
    displayTrack(igv, CUT_track)
    
    SUTs.gff3 <- read.table("./data/genome_annotations/Xu_2009_SUTs_V64.csv",
                            sep=",", as.is=TRUE)
    SUT_track <- GFF3Track("SUTs Xu 2009", SUTs.gff3,
                           url=NA_character_, indexURL=NA_character_, displayMode="EXPANDED", trackHeight=70,
                           visibilityWindow=100000)
    displayTrack(igv, SUT_track)
    
    MTX_track <- GWASTrack(glue('{input$ppi} under {input$drug} (MTX+)'), MTX_GWAS_table(), chrom.col=2, pos.col=3, pval.col=6, trackHeight=200)
    displayTrack(igv, MTX_track)
    
    noDrug_MTX_track <- GWASTrack(glue('{input$ppi} under noDrug (MTX+)'), noDrug_MTX_GWAS_table(), chrom.col=2, pos.col=3, pval.col=6, trackHeight=200)
    displayTrack(igv, noDrug_MTX_track)
    
    noDrug_noMTX_track <- GWASTrack(glue('{input$ppi} under noDrug (MTX-)'), noDrug_noMTX_GWAS_table(), chrom.col=2, pos.col=3, pval.col=6, trackHeight=200)
    displayTrack(igv, noDrug_noMTX_track)

    # noMTX_track <- GWASTrack(glue('{input$ppi} under {input$drug} (MTX-)'), noMTX_GWAS_table(), chrom.col=2, pos.col=3, pval.col=6, trackHeight=200)
    # displayTrack(igv, noMTX_track)
  
  })
}
# Create Shiny app ----

  
shinyApp(ui = ui, server = server)

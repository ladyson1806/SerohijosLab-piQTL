library(shiny)
library(plotly)
options(shiny.host = '0.0.0.0')
options(shiny.port = 3838)

# Define UI for app that draws a histogram ----

PPI_table = read.csv('./data/PPI_reference_barcodes.csv')

ui <- fluidPage(
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(width=12,
    # Input: Text box for sidplaying PPI and Drug conditions ----
      selectInput(inputId = "ppi", h4("PPI of interest"), 
                choices = PPI_table$PPI),
      selectInput(inputId = "drug", h4("Drug of interest"), 
                choices = c("Fluconazole", "5.FC", "Metformin", "Trifluoperazine"))
  ),
    
    # Main panel for displaying outputs ----
    mainPanel(width=12,
      h3("Drug + MTX"),
      br(),
      br(),
      column(
        width = 8,
        plotlyOutput(outputId = "MTX_manhattanplot")
        ),
      
      column(
        width = 4,
        plotlyOutput(outputId = "MTX_qqplot")
      ),
      
      # column(
      #   width = 3,
      #   plotlyOutput(outputId = "MTX_volcanoplot")
      # ),
      
      br(),
      br(),
      h3("noDrug + MTX"),
      column(
        width = 8,
        plotlyOutput(outputId = "noMTX_manhattanplot")
      ),
      
      column(
        width = 4,
        plotlyOutput(outputId = "noMTX_qqplot")
      ),
      
      # column(
      #   width = 3,
      #   plotlyOutput(outputId = "noMTX_volcanoplot")
      # )   
    )
  )
)

# Define server logic required to draw a histogram ----
library(glue)
library(manhattanly)

server <- function(input, output) {
  
  MTX_datasetInput <- reactive({
    PPI <- input$ppi
    DRUG <- input$drug
    infile = glue('https://raw.githubusercontent.com/ladyson1806/public_hosting/main/piQTL_mapping/formatted_for_manhattan/{PPI}_MTX_{DRUG}_avg_logratio_Fitness_minus_ref.csv')
    read.csv(infile)
  })
  
  noMTX_datasetInput <- reactive({
    PPI <- input$ppi
    DRUG <- input$drug
    infile = glue('https://raw.githubusercontent.com/ladyson1806/public_hosting/main/piQTL_mapping/formatted_for_manhattan/{PPI}_MTX_noDrug_avg_logratio_Fitness_minus_ref.csv')
    
    read.csv(infile)
  })
  
  
  output$MTX_manhattanplot <- renderPlotly({
    manhattanly(MTX_datasetInput(), chr='CHR', snp="SNP", gene="GENE", labelChr=c(1:16,"MT"), suggestiveline = -log10(1e-04), genomewideline=F, title=glue('{input$ppi} with {input$drug}' ))
  })
  
  output$MTX_qqplot <- renderPlotly({
    qqly(MTX_datasetInput(), chr='CHR', snp="SNP", gene="GENE")
  })
  
  output$MTX_volcanoplot <- renderPlotly({
    volcanoly(MTX_datasetInput(), genomewideline=-log10(1e-04))
  })
  
  output$noMTX_manhattanplot <- renderPlotly({
    manhattanly(noMTX_datasetInput(), chr='CHR', snp="SNP", gene="GENE", labelChr=c(1:16,"MT"), suggestiveline = -log10(1e-04), genomewideline=F, title=glue('{input$ppi} without {input$drug}'))
  })
  
  output$noMTX_qqplot <- renderPlotly({
    qqly(noMTX_datasetInput(), chr='CHR', snp="SNP", gene="GENE")
  })
  
  # output$noMTX_volcanoplot <- renderPlotly({
  #   volcanoly(noMTX_datasetInput(), genomewideline=-log10(1e-04))
  # })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)

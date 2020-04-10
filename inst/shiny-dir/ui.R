#' @title unCOVERApp user interface
#'
#' @description This function allows you to open graphical interface.
#' @param agree visualization.
#' @keywords app
#' @export
#' @examples
#' uncoverappLib.run()
#require(shiny)
suppressWarnings(suppressMessages(require(shiny)))
suppressWarnings(suppressMessages(require(shinyWidgets)))
suppressWarnings(suppressMessages(require(shinyBS)))
suppressWarnings(suppressMessages(require(shinyjs)))
suppressWarnings(suppressMessages(require(markdown)))
suppressWarnings(suppressMessages(require(DT)))
#options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")
suppressWarnings(suppressMessages(require(data.table)))
suppressWarnings(suppressMessages(require(dplyr)))
suppressWarnings(suppressMessages(require(Gviz)))
suppressWarnings(suppressMessages(require(Homo.sapiens)))
suppressWarnings(suppressMessages(require(OrganismDbi)))
suppressWarnings(suppressMessages(require(stringr)))
suppressWarnings(suppressMessages(require(condformat)))
suppressWarnings(suppressMessages(require(bedr)))
suppressWarnings(suppressMessages(require(openxlsx)))
suppressWarnings(suppressMessages(require(TxDb.Hsapiens.UCSC.hg19.knownGene)))
suppressWarnings(suppressMessages(require(TxDb.Hsapiens.UCSC.hg38.knownGene)))
suppressWarnings(suppressMessages(require(BSgenome.Hsapiens.UCSC.hg19)))
suppressWarnings(suppressMessages(require(BSgenome.Hsapiens.UCSC.hg38)))
suppressWarnings(suppressMessages(require(EnsDb.Hsapiens.v75)))
suppressWarnings(suppressMessages(require(EnsDb.Hsapiens.v86)))
suppressWarnings(suppressMessages(require(org.Hs.eg.db)))
#unloadNamespace("shiny")

first.file <- system.file(
  "extdata",
  "intro.md",
  package = "uncoverappLib"
)
second.file <- system.file(
  "extdata",
  "script.js",
  package = "uncoverappLib"
)
third.file <- system.file(
  "extdata",
  "CONTACTS.md",
  package = "uncoverappLib"
)

mdStyle <- "margin-left: 30px; margin-right: 30px"

intro <- function() {
  tabPanel("home",
           includeMarkdown(first.file),
           #includeMarkdown(file.path(".","intro.md")),
           style=mdStyle,
           downloadLink(outputId = "instructionsScript",
                        label="Download script for making input files",
                        style = "color:white ; background-color: #0067cd"),
           downloadLink(outputId = "dependence",
                        label="Download required dependecies for the script",
                        style = "color:white ; background-color: #0067cd"),
           downloadLink(outputId = "configuration_file",
                        label= "Download script for configration bash scritp",
                        style = "color:white ; background-color: #0067cd"),
           br(),
           br(),
           br()

  )
}
myHome <- function() {
  tabPanel("Coverage Analysis",
           h1(strong("Interactive web-application to visualize and
                     annotate low-coverage positions in clinical sequencing")),
           #titlePanel("Coverage sequencing Explorer"),
           helpText(em("Note:Select input options",
                       span("Upload your input bed.gz
                            file with columns: chromosome, start, end, coverage
                            by sample", style = "color:blue"))),
           shinyjs::useShinyjs(),
           includeScript(second.file),
           sidebarPanel(
             selectInput("UCSC_Genome",
                         label = "Reference Genome",
                         choices = c("hg19",
                                     "hg38"),
                         selected = "UCSC genome"),
             hr(),

             textInput(inputId = "Gene_name",
                       label= "Gene name"),
             #actionButton("button1",label= "Apply"),
             bsButton("button1",label= "Apply",  icon = icon("power-off"),
                      style = "success", size = "extra-small"),
             shinyjs::hidden(p(id = "text1", "Running.....")),
             #actionButton("remove",label= "Refresh"),
             bsButton("remove",label= "Refresh",  icon = icon("power-off"),
                      style = "success", size = "extra-small"),
             helpText(em("write gene name corrisponding coordinate
                         positions and action button apply")),
             hr(),

             pickerInput("Chromosome",
                         label = "Chromosome",
                         choices = c("chr1", "chr2","chr3", "chr4","chr5",
                                     "chr6","chr7","chr8","chr9","chr10",
                                     "chr11","chr12","chr13","chr14","chr15",
                                     "chr16","chr17",
                                     "chr18","chr19", "chr20", "chr21",
                                     "chr22", "chrX", "chrY","chrM",
                                     names("file1")),
                         options = list(`actions-box` = TRUE),
                         multiple =FALSE),
             hr(),
             pickerInput("coverage_co",
                         label = "Coverage threshold",
                         choices = c(1:1000, "all",names("file1")),
                         options = list(`actions-box` = TRUE), multiple =FALSE),
             helpText(em("Select minimum value as coverage threshold")),
             hr(),
             textInput(inputId = "Sample",
                       label= "Sample"),

             helpText(em("Select sample for coverage analysis.
                         Example:sample_1")),

             hr(),
             splitLayout(cellWidths = c("30%", "70%"),
                         textInput(inputId = "transcript_id",
                                   label= "Transcript number"),
                         textInput(inputId = "id_t",
                                   label= "Transcript ID")),

             helpText(em("Retrieve your favourite transcript number from UCSC exons")),


             hr(),
             pickerInput("exon_number",
                         label= "exon number",
                         choices = c(1:150),
                         options = list(`actions-box` = TRUE), multiple =FALSE),
             #actionButton("button5",label= "Make exon"),
             bsButton("button5",label= "Make exon",
                      icon = icon("power-off"), style = "success",
                      size = "extra-small"),
             shinyjs::hidden(p(id = "text1", "Running.....")),
             #actionButton("remove5",label= "Refresh"),
             bsButton("remove5",label= "Refresh",
                      icon = icon("power-off"), style = "default",
                      size = "extra-small"),
             helpText(em("zooming one exon")),
             hr(),
             hr(),

             textInput(inputId = "Start_genomicPosition",
                       label = "START genomic position"),

             textInput(inputId = "end_genomicPosition",
                       label = "END genomic position"),
             helpText(em("change genomic intervall for zooming")),


             hr(),
             hr(),
             textInput(inputId = "query_Database",
                       label= "Region coordinates"),


             helpText(em("write to expland dbNSFP-annotated genomic positions.
                         For example 2:166845670-166930180")),
             hr(),

             hr(),
             downloadButton("downloadData", "Download", class = "btn-primary",
                            style='padding:4px; font-size:120%'),
             hr(),
             hr(),
             shiny::tags$button(
               id = 'close',
               type = "button",
               class = "btn action-button",
               style='color: white; background-color: #dd4b39;
               padding:4px; font-size:120%',
               onclick = "setTimeout(function(){window.close();},500);",
               # close browser
               "Close App")),
           mainPanel(
             fileInput(inputId = "file1",
                       label = "Select input file",
                       accept = c("text/csv",
                                  ".bedGraph",
                                  ".bedGraph.gz",
                                  ".zip",
                                  ".gz",
                                  ".bed",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             checkboxInput("header", "Header", TRUE),
             shinyBS::bsButton("example_data",
                               label= "load example dataset",
                               icon = icon("file"), style = "info",
                               size = "extra-small"),
             helpText(em("load example dataset with
                         base coverage counts of POLG gene")),
             hr(),

             tabsetPanel(
               tabPanel("bed file", DT::dataTableOutput("text")),
               tabPanel("UCSC gene", DT::dataTableOutput("ccg")),
               tabPanel("UCSC exons",
                        helpText(em("Extract protein coding positions from
                                    UCSC")), DT::dataTableOutput("exon_pos")),
               tabPanel("Low-coverage positions",
                        DT::dataTableOutput("text_cv")),
               tabPanel("Gene coverage", plotOutput("all_gene"),
                        DT::dataTableOutput('df.l')),
               tabPanel("Exon coverage", helpText(em("This is a Gviz function
               and it plots exon with ideogram, genome coordinates,
               coverage information, Ensembl and UCSC gene annotation.The
              annotation for the databases are directly fetched from Ensemb and
                      all tracks will be plotted in a 3' -> 5' direction.")),
                        plotOutput("ens",
                                   dblclick = "plot1_dblclick",
                                   brush = brushOpts( id = "plot1_brush",
                                                      resetOnNew = TRUE)),
                        DT::dataTableOutput("tabExon")),
               tabPanel("Zoom to sequence",
                        helpText(em("choose a few genomic intervals")),
                        plotOutput("sequence")),
               tabPanel("Annotations on low-coverage positions",
                        helpText(em("dbSNP-annotation collects all
                                    consequences found in VEP-defined
                                    canonical transcripts")),
                        condformatOutput("uncover_position"))
             ),
             hr(),

             fluidRow(
               column(12,DT::dataTableOutput("x4"))
             )
           )
  )
}
myTab1 <- function() {
  tabPanel("Calculate AF by allele frequency app",

           # Application title
           titlePanel("Maximum credible population allele frequency"),

           ##### Bootstrap method for page costruction
           fluidPage(
             fluidRow(
               ##### Sidebar
               column(8,wellPanel(radioButtons("inh",
                                               "Inheritance:",
                                               choices = list("monoallelic",
                                                              "biallelic"),
                                               selected = "monoallelic"),

                                  numericInput("prev","Prevalence = 1 in ...
                                               (people)",
                                               min = 1,max = 1e8,value = 500),
                                  options = NULL),
                      br(),
                      sliderInput("hetA","Allelic heterogeneity:",
                                  min = 0, max = 1,value = 0.1),
                      sliderInput("hetG",
                                  "Genetic heterogeneity:",
                                  min = 0, max = 1,value = 1),
                      br(),
                      sliderInput("pen", "Penetrance:",
                                  min = 0, max = 1, value = 0.5))),
             br(),
             column(8,
                    h3("Maximum credible population AF:"),
                    h2(textOutput("maxAF"),align="center",style = "color:blue")),
             column(8,
                    h3("Uncorver position",
                       helpText(em("Low-coverage positions excluding sites
                                   annotated as variants with AF> maxAF
                                   (default maxAF value: 5%)"),align="center",
                                style="color:blue"),
                       style = "font-size: 100%; width: 100%",
                       condformatOutput("uncoverPosition"))),
             br(),
             br(),
             downloadButton("download_maxAF", "Download_maxAF",
                            class = "btn-primary",
                            style='padding:4px; font-size:80%',
                            helpText("download low coverage
                                     position dbSNFP-annotation filtered by
                                     maximum allele frequency",
                                     class = "btn-primary",
                                     style='padding:4px; font-size:60%'))
             #)
           ))}


myTab2 <- function() {
  tabPanel("Binomial distribution",
           titlePanel("Binomial distribution "),
           fluidRow(
             column(4,(numericInput("p",
                                    "Allele Fraction",
                                    min = 0,
                                    max = 1,
                                    value = 0.00001)),
                    hr(),
                    numericInput("num_all",
                                 "Variant reads",
                                 min=0,
                                 max=100000,
                                 value=10),
                    hr(),

                    textInput(inputId = "start_gp",
                              label = "START genomic position"),

                    textInput(inputId = "end_gp",
                              label = "END genomic position"),
                    helpText(em("Specify start and end coordinates
                                for your genomic region of interest"))),


             column(4,
                    h2("consideration:"),
                    h3(htmlOutput("ci"))),
             column(4,
                    h2("Binomial Distribution", plotOutput("bd"))),
             column(10, h2("Cumulative distribution function"),
                    h3(plotOutput("pbinom"))))

  )}


myabout <- function() {
  tabPanel("contacts",
           includeMarkdown(third.file),
           style=mdStyle
  )
}



ui <- shinyUI(navbarPage("uncoverApp",
                         intro(),
                         myHome(),
                         myTab1(),
                         myTab2(),
                         myabout()
                         # preprocess()
))


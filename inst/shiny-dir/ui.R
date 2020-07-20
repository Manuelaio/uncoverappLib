#' @title unCOVERApp user interface
#'
#' @description This function allows you to open graphical interface.
#' @param agree visualization.
#' @keywords app
#' @export
#' @examples
#' uncoverappLib.run()

require(shiny)
require(shinyWidgets)
require(shinyBS)
require(shinyjs)
require(markdown)
require(DT)
#options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")
#require(data.table)
#require(dplyr)
require(Gviz)
require(Homo.sapiens)
require(OrganismDbi)
require(stringr)
require(condformat)
require(shinyjs)
require(shinycssloaders)
require(bedr)
require(Rsamtools)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(BSgenome.Hsapiens.UCSC.hg19)
require(BSgenome.Hsapiens.UCSC.hg38)
require(EnsDb.Hsapiens.v75)
require(EnsDb.Hsapiens.v86)
require(org.Hs.eg.db)
options("browser" = "xdg-open")

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

intro_processing <- system.file(
  "extdata",
  "prep_input.md",
  package = "uncoverappLib"
)

mdStyle <- "margin-left: 30px; margin-right: 30px"

intro <- function() {
  tabPanel("home",
           shiny::includeMarkdown(first.file),
           #includeMarkdown(file.path(".","intro.md")),
           style=mdStyle

  )
}

preprocess <- function() {
  shiny::tabPanel("Processing",
                  shiny::includeMarkdown(intro_processing),
                  #shiny::includeMarkdown(file.path(".", "prep_input.md")),
                  h1(strong("Prepare your input file")),
                  fluidPage(
                    sidebarLayout(
                      sidebarPanel(
                        shiny::selectInput("Genome",
                                           label = "Reference Genome",
                                           choices = c("hg19",
                                                       "hg38"),
                                           selected = "UCSC genome"),
                        hr(),
                        shinyWidgets::pickerInput("notation",
                                                  label = "Chromosome Notation",
                                                  choices = c("chr", "number"),
                                                  options = list(`actions-box` = TRUE),
                                                  multiple =FALSE),

                        shinyWidgets::pickerInput("MAPQ",
                                                  label = "minum Mapping Quality (MAPQ)",
                                                  choices = c(1:1000),
                                                  options = list(`actions-box` = TRUE), multiple =FALSE),

                        shinyWidgets::pickerInput("base_qual",
                                                  label = "minimum Base Quality",
                                                  choices = c(1:1000),
                                                  options = list(`actions-box` = TRUE), multiple =FALSE),
                        hr(),

                        fileInput(inputId = "gene1",
                                  label = "Load a txt file with gene(s) list: one gene
                              in one row",
                                  accept = c("text/csv",
                                             ".zip",
                                             ".gz",
                                             ".bed",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),

                        fileInput(inputId = "bam_list",
                                  label = "Load a bam.list file: one bam paths
                                in one row",
                                  accept = c("text/csv",
                                             ".zip",
                                             ".gz",
                                             ".bed",
                                             "text/comma-separated-values,text/plain",
                                             ".list"))),
                      shiny::mainPanel(
                        tabsetPanel(
                          tabPanel(title= "input for uncoverapp",
                                   shinycssloaders::withSpinner(
                                     dataTableOutput("input1"))
                          )
                        ) )
                    )))
}

myHome <- function() {
  tabPanel("Coverage Analysis",
           h1(strong("Interactive web-application to visualize and
                     annotate low-coverage positions in clinical sequencing")),
           #titlePanel("Coverage sequencing Explorer"),
           helpText(em("Note:Select input options",
                       span("Upload your input bed
                            file with columns: chromosome, start, end, coverage
                            by sample and nucleotide count", style = "color:blue"))),
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
             shinyBS::bsButton("button1",label= "Apply",  icon = icon("power-off"),
                               style = "success", size = "extra-small"),
             shinyjs::hidden(p(id = "text1", "Running.....")),
             #actionButton("remove",label= "Refresh"),
             shinyBS::bsButton("remove",label= "Refresh",  icon = icon("power-off"),
                               style = "success", size = "extra-small"),
             helpText(em("write gene name and push apply button")),
             hr(),

             textInput(inputId ="Chromosome",
                        label = "Chromosome"),

             #shinyWidgets::pickerInput("Chromosome",
              #                         label = "Chromosome",
               #                        choices = c("chr1", "chr2","chr3", "chr4","chr5",
                #                                   "chr6","chr7","chr8","chr9","chr10",
                 #                                  "chr11","chr12","chr13","chr14","chr15",
                  #                                 "chr16","chr17",
                   #                                "chr18","chr19", "chr20", "chr21",
                    #                               "chr22", "chrX", "chrY","chrM",
                     #                              names("file1")),
                      #                 options = list(`actions-box` = TRUE),
                       #                multiple =FALSE),
             hr(),
             shinyWidgets::pickerInput("coverage_co",
                                       label = "Coverage threshold",
                                       choices = c(1:1000, "all",names("file1")),
                                       options = list(`actions-box` = TRUE), multiple =FALSE),
             helpText(em("Select minimum value as coverage threshold")),
             hr(),
             textInput(inputId = "Sample",
                       label= "Sample"),

             helpText(em("Select number of sample for coverage analysis.
                         Example:1")),

             hr(),
             splitLayout(cellWidths = c("30%", "70%"),
                         textInput(inputId = "transcript_id",
                                   label= "Transcript number"),
                         textInput(inputId = "id_t",
                                   label= "Transcript ID")),

             helpText(em("Retrieve your favourite transcript number from UCSC exons")),


             hr(),
             shinyWidgets::pickerInput("exon_number",
                                       label= "exon number",
                                       choices = c(1:150),
                                       options = list(`actions-box` = TRUE), multiple =FALSE),
             #actionButton("button5",label= "Make exon"),
             shinyBS::bsButton("button5",label= "Make exon",
                               icon = icon("power-off"), style = "success",
                               size = "extra-small"),
             shinyjs::hidden(p(id = "text1", "Running.....")),
             #actionButton("remove5",label= "Refresh"),
             shinyBS::bsButton("remove5",label= "Refresh",
                               icon = icon("power-off"), style = "default",
                               size = "extra-small"),
             helpText(em("zooming one exon")),
             hr(),
             hr(),

             textInput(inputId = "Start_genomicPosition",
                       label = "START genomic position"),

             textInput(inputId = "end_genomicPosition",
                       label = "END genomic position"),
             helpText(em("change genomic interval for zooming")),


             hr(),
             hr(),
             textInput(inputId = "query_Database",
                       label= "Region coordinates"),


             helpText(em("write to expand dbNSFP-annotated genomic positions.
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
             shinyBS::bsButton("pileup",
                               label= "load input file",
                               icon = icon("file"), style = "info",
                               size = "default"),
             #shinyBS::bsButton("example_data",
             #                 label= "load example dataset",
             #                icon = icon("file"), style = "info",
             #               size = "extra-small"),
             #helpText(em("load example dataset with
             #           base coverage counts of POLG gene")),
             hr(),

             tabsetPanel(
               tabPanel("bed file", shinycssloaders::withSpinner(
                 DT::dataTableOutput("text"))),
               tabPanel("UCSC gene",
                        DT::dataTableOutput("ccg")),
               tabPanel("UCSC exons",
                        helpText(em("Extract protein coding positions from
                                    UCSC")), DT::dataTableOutput("exon_pos")),
               tabPanel("Low-coverage positions",
                        DT::dataTableOutput("text_cv")),
               tabPanel("Gene coverage", shinycssloaders::withSpinner(
                 plotOutput("all_gene")),
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
                        shinycssloaders::withSpinner(
                          condformat::condformatOutput("uncover_position"))),
               id = "tabSet"
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
                    h3("Uncover position",
                       helpText(em("Low-coverage positions excluding sites
                                   annotated as variants with AF> maxAF
                                   (default maxAF value: 5%)"),align="center",
                                style="color:blue"),
                       style = "font-size: 100%; width: 100%",
                       shinycssloaders::withSpinner(
                       condformat::condformatOutput("uncoverPosition")))),
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
                                    value = 0.05)),
                    helpText(em("the expected fraction of variant reads,
                    probability of success",
                    align="center",
                    style="color:gray")),
                    hr(),
                    numericInput("num_all",
                                 "Variant reads",
                                 min=0,
                                 max=100000,
                                 value=10),

                    helpText(em("the minimum number of variant reads required
                                by the user to support variant calling,
                                (number of successes)"),align="center",
                             style="color:gray"),
                    hr(),

                    textInput(inputId = "start_gp",
                              label = "START genomic position"),

                    textInput(inputId = "end_gp",
                              label = "END genomic position"),
                    helpText(em("Specify start and end coordinates
                                for your genomic region of interest gene"))),


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
                         preprocess(),
                         myHome(),
                         myTab1(),
                         myTab2(),
                         myabout()
))


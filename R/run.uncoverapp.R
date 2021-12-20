#' @title run.uncoverapp
#'
#' @description This function launches \code{unCOVERApp}, a
#' \code{Shiny} application for clinical assessment of sequence coverage.
#' Setting where uncoverapp will be launched with following \code{where} option:
#' `"browser`" in user default browser,  `"panel`"  RStudio viewer and `"window`"
#' in a new Rstudio window.
#'
#' @author Emanuela Iovino
#'
#'
#' @examples
#' \dontrun{
#' file.name='../path/sorted.bed.gz'
#' tbi='.../path/sorted.bed.gz.tbi'
#' run.uncoverapp(where="window")
#' }
#'
#' ## Only run this example in interactive R sessions
#'
#' if (interactive()) {
#' run.uncoverapp(where="window")
#' }
#'
#' #After running `run.uncoverapp(where="window")` the shiny app appears in your choosen location.
#'
#'
#' @return
#'  This return a Shiny App. The is no value
#
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @rawNamespace import(Gviz, except = tags)
#' @import shinyWidgets
#' @import shinyBS
#' @importFrom shinyjs useShinyjs hidden enable
#' @importFrom shinycssloaders withSpinner
#' @import markdown
#' @importFrom DT renderDataTable dataTableOutput
#' @import Homo.sapiens
#' @import openxlsx
#' @import condformat
#' @import stringr
#' @import processx
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Hsapiens.v75
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import OrganismDbi
#' @import tools
#' @importFrom  Rsamtools ScanBamParam PileupParam pileup
#' @importFrom rlist list.append
#' @importFrom GenomicRanges makeGRangesFromDataFrame

#'
#' @export
#' @param where accept `"browser`" , `"panel`"  or `"window`". The option sets
#' where uncoverapp will be launched

run.uncoverapp <- function(where=c("browser", "panel", "window")){
  message("uncoverapp must be called in RStudio")
  options(shiny.launch.browser = switch(
    match.arg(where, c("browser", "pane", "window")),
    browser = get(".rs.invokeShinyWindowExternal", "tools:rstudio"),
    panel = get(".rs.invokeShinyPaneViewer", "tools:rstudio"),
    window = get(".rs.invokeShinyWindowViewer", "tools:rstudio")
  ))
  appDir <- system.file("shiny-dir", package = "uncoverappLib")
  if (appDir == "") {
    stop("Try reinstalling uncoverappLib .", call = FALSE)
  }

  shiny::runApp(appDir)
}

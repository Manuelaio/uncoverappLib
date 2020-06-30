#' @title run.uncoverapp
#'
#' @description This function launches \code{unCOVERApp}, a
#' \code{Shiny} application for clinical assessment of sequence coverage.
#' @author Emanuela Iovino
#'
#'
#' @examples
#' \donttest{
#' file.name='../path/sorted.bed.gz'
#' tbi='.../path/sorted.bed.gz.tbi'
#' run.uncoverapp()
#' }
#'
#' ## Only run this example in interactive R sessions
#'
#' if (interactive()) {
#' run.uncoverapp()
#' }
#'
#' #After running `run.uncoverapp()` the shiny app appears in your browser
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
#' @importFrom  Rsamtools ScanBamParam PileupParam pileup
#' @importFrom rlist list.append
#' @importFrom GenomicRanges makeGRangesFromDataFrame

#'
#' @export
#'

run.uncoverapp <- function() {

  appDir <- system.file("shiny-dir", package = "uncoverappLib")
  if (appDir == "") {
    stop("Try reinstalling uncoverappLib .", call = FALSE)
  }

  shiny::runApp(appDir)
}

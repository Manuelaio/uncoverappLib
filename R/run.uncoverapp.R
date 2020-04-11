#' @title run.uncoverapp
#'
#' @description This function launches \code{unCOVERApp}, a \code{Shiny} application for clinical assessment of sequence coverage.
#' @author Emanuela Iovino
#'
#'
#' @examples
#' \dontrun{
#' run.uncoverapp()
#' to run unCOVERApp use run.uncoverapp() function.
#' }
#'
#'
#' #After running `run.uncoverapp()` the shiny app appears in your browser
#' @return
#'  This return a Shiny App. The is no value
#
#' @import shiny
#' @import shinyWidgets
#' @import shinyBS
#' @importFrom shinyjs useShinyjs hidden enable
#' @import markdown
#' @import dplyr
#' @import Homo.sapiens
#' @import stringr
#' @import condformat
#' @import bedr
#' @import openxlsx
#' @import readxl
#' @import writexl
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import EnsDb.Hsapiens.v75
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Hsapiens.v86
#' @import org.Hs.eg.db
#' @import data.table
#' @import OrganismDbi
#' @importFrom Gviz SequenceTrack DataTrack IdeogramTrack
#' @importFrom Gviz GeneRegionTrack OverlayTrack GenomeAxisTrack
#' @importFrom Gviz plotTracks
#' @importFrom grDevices extendrange
#' @importFrom DT renderDataTable dataTableOutput
#'
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

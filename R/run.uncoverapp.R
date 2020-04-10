#' @title uncoverappLib.run
#'
#' @description This function launches \code{unCOVERApp}, a \code{Shiny} application for clinical assessment of sequence coverage.
#' @author Emanuela Iovino
#'
#'
#' @keywords uncoverappLib
#'
#' @examples
#' #to run unCOVERApp use uncoverappLib.run() function.
#'
#' #run.uncoverapp()
#'
#' #After running `uncoverappLib.run()`, users use *unCOVERApp*
#' @return
#'  This return `Shiny` App. The is no value
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
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import EnsDb.Hsapiens.v75
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Hsapiens.v86
#' @import org.Hs.eg.db
#'
#'
#'
#' @export
#'
#'
#

run.uncoverapp <- function() {

  appDir <- system.file("shiny-dir", package = "uncoverappLib")
  if (appDir == "") {
    stop("Try reinstalling uncoverappLib .", call = FALSE)
  }

  shiny::runApp(appDir)
}

#' @title run.uncoverapp
#'
#' @description This function launches \code{unCOVERApp}, a
#' \code{Shiny} application for clinical assessment of sequence coverage.
#' Setting where uncoverapp will be launched with following \code{where} option:
#' `"browser`" in user default browser,  `"viewer`"  RStudio viewer and `"window`"
#' in a new Rstudio window.
#'
#' @author Emanuela Iovino
#'
#'
#' @examples
#' \dontrun{
#' file.name='../path/sorted.bed.gz'
#' tbi='.../path/sorted.bed.gz.tbi'
#' app()
#' }
#'
#' ## Only run this example in interactive R sessions
#'
#' if (interactive()) {
#' app()
#' }
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
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @import OrganismDbi
#' @importFrom  Rsamtools ScanBamParam PileupParam pileup
#' @importFrom rlist list.append
#' @importFrom GenomicRanges makeGRangesFromDataFrame

#'
#' @export


uncoverAPP<- function(){
  appDir <- system.file("shiny-dir", package = "uncoverappLib")
  if (appDir == "") {
    stop("Try reinstalling uncoverappLib .", call = FALSE)
  }
  shiny::runApp(appDir)
}


#'  Location for uncoverapp in RStudio enviroment
#'
#' This function controls the `shiny.launch.browser` option to launch uncoverapp
#' in an external `browser`, the RStudio viewer `"viewer"`, or a new `"window"` in
#' RStudio.
#' @param where accept `"browser`" , `"viewer`"  or `"window`". The option sets
#' where uncoverapp will be launched. Using NULL, uncoverapp will use default
#' After running `run.uncoverapp(where="window")` the shiny app appears
#' in your chosen location.
#'
#' @return
#' This return a  Shiny App. The is no value
#'
#' @examples
#' ## Only run this example in interactive R sessions
#'
#' if (interactive()) {
#' run.uncoverapp(where="window")
#' }
#'
#'
#' @export


run.uncoverapp <- function(where=c("browser", "viewer", "window")){
  if(missing(where)){
    message("uncoverapp is launching with your default option")}else{message(
      "where function must be called in RStudio"
    )}
  if(missing(where)){
    uncoverAPP()
  }else{
    options(shiny.launch.browser = switch(
      match.arg(where, c("browser", "viewer", "window")),
      browser = get(".rs.invokeShinyWindowExternal", "tools:rstudio"),
      viewer = get(".rs.invokeShinyPaneViewer", "tools:rstudio"),
      window = get(".rs.invokeShinyWindowViewer", "tools:rstudio")
    ))
    uncoverAPP()}
}



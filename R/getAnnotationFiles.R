# functions to download annotation file from Zenodo for uncoverappLib

#' wrapper function for getting BiocFileCache associated with uncoverapp package
#'
#' @return BiocFileCache object associated with uncoverappLib
#' @import BiocFileCache
#' @importFrom rappdirs user_cache_dir
.get_cache <- function() {
  cache <- rappdirs::user_cache_dir(appname = "uncoverappLib")
  BiocFileCache::BiocFileCache(cache,ask=FALSE)
}

#' download and rename sorted.bed.gz and sorted.bed.gz.tbi files for annotation of
#' low-coverage positions.
#'
#' @param verbose (logical) print messages
#' @examples getAnnotationFiles()
#' @return (char) Path to local cached file
#' or initial download is required
#' @export
#'
#'
getAnnotationFiles <- function(verbose = FALSE) {
  fileURL <- "https://zenodo.org/record/3747448/files/sorted.bed.gz"
  fileURL2 <- "https://zenodo.org/record/3747448/files/sorted.bed.gz.tbi"
  bfc <- .get_cache()
  rid <- bfcquery(bfc, "sorted.bed.gz$", "rname")$rid
  rid2 <- bfcquery(bfc, "sorted.bed.gz.tbi", "rname")$rid
  bfcrpath(bfc, rids = rid)

  if (!length(rid)) {
    if( verbose )
      message( "Downloading GENE file" )
    rid <- names(bfcadd(bfc, "sorted.bed.gz", fileURL ))

  }

  if (!isFALSE(bfcneedsupdate(bfc, rid)))
    bfcdownload(bfc, rid)


  bfcrpath(bfc, rids = rid2)
  if (!length(rid2)) {
    if( verbose )
      message( "Downloading tbi file" )
    rid2 <- names(bfcadd(bfc, "sorted.bed.gz.tbi", fileURL2 ))
  }

  if (!isFALSE(bfcneedsupdate(bfc, rid2)))
    bfcdownload(bfc, rid2)


  #bfcrpath(bfc, rids = c(rid,rid2))

  bfcrpath(bfc, rids = c(rid,rid2))
  path_out=base::file.path(bfcrpath(bfc, rids= rid))

  rename1=base::gsub('[[:digit:]].*', 'sorted.bed.gz', path_out)

  path_tbi= base::file.path(bfcrpath(bfc, rids= rid2))
  rename2=base::gsub('[[:digit:]].*', 'sorted.bed.gz.tbi', path_tbi)

  base::file.copy(from = path_out, to = rename1)
  base::file.copy(from = path_tbi, to = rename2)
  return(print(c(rename1, rename2)))

}




# functions to download annotation file from Zenodo for uncoverappLib

#' wrapper function for getting BiocFileCache associated with uncoverapp package
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
#'
#' @param verbose (logical) print messages
#' @param vignette (logical) download example annotation-file in vignette mode
#' @examples
#' getAnnotationFiles(verbose = TRUE, vignette= TRUE)
#'
#' @return (char) Path to local cached file
#' or initial download is required
#' @export
#'
#'
getAnnotationFiles <- function(verbose = FALSE, vignette= FALSE) {
  if (vignette == FALSE){
  fileURL <- "https://zenodo.org/record/3747448/files/sorted.bed.gz"
  fileURL2 <- "https://zenodo.org/record/3747448/files/sorted.bed.gz.tbi"
  bfc <- .get_cache()
  rid <- bfcquery(bfc, "sorted.bed.gz$", "rname")$rid

  rid2 <- bfcquery(bfc, "sorted.bed.gz.tbi", "rname")$rid
  bfcrpath(bfc, rids = rid)

  if (!length(rid)) {
    if( verbose )
      message( "Downloading annotations file, please wait few minutes if you
      launch getAnnotationFiles() for the first time" )
    rid <- names(bfcadd(bfc, "sorted.bed.gz", fileURL ))
  }
  if(length(rid)){
    message("your file already exists in cache")
  }

  #if (!isFALSE(bfcneedsupdate(bfc, rid)))
   # bfcdownload(bfc, rid)


  bfcrpath(bfc, rids = rid2)
  if (!length(rid2)) {
    if( verbose )
      message( "Downloading tbi file" )
    rid2 <- names(bfcadd(bfc, "sorted.bed.gz.tbi", fileURL2 ))
  }

  #if (!isFALSE(bfcneedsupdate(bfc, rid2)))
   # bfcdownload(bfc, rid2)

  bfcrpath(bfc, rids = c(rid,rid2))

  path_out=base::file.path(bfcrpath(bfc, rids= rid))

  rename1=base::gsub('[[:digit:]].*', 'sorted.bed.gz', path_out)

  path_tbi= base::file.path(bfcrpath(bfc, rids= rid2))
  rename2=base::gsub('[[:digit:]].*', 'sorted.bed.gz.tbi', path_tbi)

  if(!file.exists(rename1[1])){
    if( verbose )
      message( "Rename file in cache, please wait few minutes" )
  base::file.copy(from = path_out, to = rename1)

  }

  if(!file.exists(rename2[1])){
  base::file.copy(from = path_tbi, to = rename2)
    }
  return(print(c(rename1, rename2)))
  }else{
    example_bfc <- .get_cache()

    exampleDATA<-"https://zenodo.org/record/3909001/files/POLG.bed.gz"
    exampleTBI<-"https://zenodo.org/record/3909001/files/POLG.bed.gz.tbi"
    example_rid<- bfcquery(example_bfc, "POLG.bed.gz$", "rname")$rid
    example_rid2<- bfcquery(example_bfc, "POLG.bed.gz.tbi$", "rname")$rid
    bfcrpath(example_bfc, rids = example_rid)
    if (!length(example_rid)) {
      if( verbose )
        message( "Downloading annotations file, please wait few minutes if you
      launch getAnnotationFiles() for the first time" )
      example_rid <- names(bfcadd(example_bfc, "POLG.bed.gz", exampleDATA ))
    }
    if(length(exampleDATA)){
      message("your file already exists in cache")
    }
    bfcrpath(example_bfc, rids = example_rid2)
    if (!length(example_rid2)) {
      if( verbose )
        message( "Downloading tbi file" )
      example_rid2 <- names(bfcadd(example_bfc, "POLG.bed.gz.tbi", exampleTBI ))
    }

    bfcrpath(example_bfc, rids = c(example_rid,example_rid2))

    example_path_out=base::file.path(bfcrpath(example_bfc, rids= example_rid))

    example_rename1=base::gsub('[[:digit:]].*', 'POLG.bed.gz', example_path_out)

    example_path_tbi= base::file.path(bfcrpath(example_bfc, rids= example_rid2))
    example_rename2=base::gsub('[[:digit:]].*', 'POLG.bed.gz.tbi', example_path_tbi)

    if(!file.exists(example_rename1[1])){
      if( verbose )
        message( "Rename file in cache, please wait few minutes" )
      base::file.copy(from = example_path_out, to = example_rename1)

    }

    if(!file.exists(example_rename2[1])){
      base::file.copy(from = example_path_tbi, to = example_rename2)
    }
    return(print(c(example_rename1, example_rename2)))
  }

}




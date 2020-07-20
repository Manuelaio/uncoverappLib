#'Build input file
#'
#' @description
#' Function to build input file for unCOVERAPP when the number
#' of genes to analyze is > 50.
#'
#'
#' @return A file.bed containing tab-separated specifications of
#' genomic coordinates (chromosome, start position, end position),
#' the coverage value, and the reference:alternate allele counts for each position.
#'
#' @param geneList a text file, named with .txt extension,
#' containing HGNC official gene name(s) one per row.
#' @param genome (char) reference genome, hg19 or hg38
#' @param type_bam (char) chromosome notation of their BAM file(s). Use "number"
#' or "chr".  In the BAM file, the number option refers to 1, 2, ..., X,.M
#' chromosome notation, while the chr option refers to chr1, chr2, ... chrX,
#' chrM chromosome notation.
#' @param bamList a text file, named with .list extension,
#' containing the absolute paths to BAM file(s) one per row.
#' @param outDir (char) directory where pileup output will be stored
#' @param MAPQ.min (integer) minimum MAPQ value for an alignment
#' to be included in output file.
#' @param base.quality (integer) minimum QUAL value for each
#' nucleotide in an alignment.
#' @export
#' @import rlist
#' @import OrganismDbi
#' @import Rsamtools
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom utils write.table
#' @import GenomicRanges
#' @examples
#' gene.list<- system.file("extdata", "mygene.txt", package = "uncoverappLib")
#'
#' bam_example <- system.file("extdata", "example_POLG.bam",
#' package = "uncoverappLib")
#' cat(bam_example, file = "bam.list", sep = "\n")
#' temp_dir=tempdir()
#' buildInput(geneList= gene.list, genome= "hg19", type_bam= "chr",
#' bamList= "bam.list", outDir= temp_dir)
#'


buildInput<- function(geneList,genome,type_bam,bamList,outDir,
                      MAPQ.min=1,base.quality=1) {


    ##### Check input
    if (missing(geneList)) stop("geneList must be supplied.\n")
    if (missing(bamList)) stop("bamList must be supplied.\n")
    if (MAPQ.min <= 0 )
    stop("MAPQ.min must be greater than 0")
    if (base.quality <= 0 )
    stop("base.quality must be greater than 0")

    ##load packages

    ###read gene(s) list and retrieve genomic coordinates in grange format

    message("Reading gene name started at:")
    message(Sys.time())
    gene.List=scan(geneList, character(), quote = "")
    my_gene_name=OrganismDbi::select(org.Hs.eg.db, key= gene.List,
                                     columns="ENTREZID", keytype="ALIAS")

    ###check control: stop when gene names in list are incorrect

    for (i in my_gene_name$ENTREZID){
     if (is.na(i)){
       s=print(subset(my_gene_name$ALIAS,is.na(my_gene_name$ENTREZID)))
       utils::write.table(s, file='./preprocessing_log.txt', quote= FALSE,
                  row.names = FALSE, col.names = FALSE)
      stop('ERROR: unrecognized GENE NAME. Please choose a
         correct HGNC official gene names.')}
    else {return= NULL}
  }
  ID=my_gene_name$ENTREZID
  if (genome == "hg19"){
    all_gene= TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene}else{
      all_gene= TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene}

  b=c()
  for (i in ID){
    if( ! is.element(i, keys(all_gene, "GENEID")))
      b[i]= i
  }
  no_entrID= subset(my_gene_name, my_gene_name$ENTREZID %in% b)

  ###chech controls : remove genes name not recognize from
  #TxDb.Hsapiens.UCSC.hg19.knownGene
  message("Control HGNC gene name started at:")
  message(Sys.time())

  if (nrow(no_entrID)!=0){
    utils::write.table(no_entrID, file='./preprocessing_log1.txt',
                quote= FALSE, row.names = FALSE, col.names = FALSE)
    stop('ERROR: unrecognized GENE NAME. Please remove genes
       stored in preprocessing_log1.txt')}


  pre= data.frame()
  for (i in ID){
    txid <- OrganismDbi::select(all_gene, i, "TXNAME", "GENEID")[["TXNAME"]]
    cds <-OrganismDbi::exonsBy(all_gene, by="tx", use.names=TRUE)
    coordinate= as.data.frame(cds[names(cds) %in% txid])
    pre= rbind(coordinate,pre)
    pre$start= pre$start -10
    pre$end= pre$end +10
  }


  for_bed=data.frame()
  #for (row in 1:nrow(pre)){
  for (row in seq_len(nrow(pre))){
    sel_cc= data.frame(as.character(pre$seqnames[row]),
                       pre$start[row], pre$end[row], stringsAsFactors = FALSE)
    for_bed= rbind(for_bed, sel_cc)
  }
  colnames(for_bed)=c('chr', 'start', 'end')

  for_bed= unique(for_bed)



  if(type_bam == "number"){
    for_bed$chr= gsub("^.{0,3}", "", for_bed$chr, perl =TRUE)}



  for_grange=GenomicRanges::makeGRangesFromDataFrame(for_bed, keep.extra.columns = TRUE)

  ##set parameters for bam

  param <- Rsamtools::ScanBamParam(which= for_grange)

  p_param <- Rsamtools::PileupParam(distinguish_nucleotides=TRUE,
                         distinguish_strands=FALSE,
                         min_mapq=MAPQ.min,
                         min_base_quality=base.quality,
                         min_nucleotide_depth=0,
                         max_depth=15000)
  list= scan(bamList, character(), quote= "")
  df= list()
  for (i in list){
    pileup_df= Rsamtools::pileup(list, scanBamParam=param, pileupParam=p_param)

    df=rlist::list.append(df, pileup_df)
  }
  #remove which_label column
  lst1 <- lapply(df, function(x) transform(x[,-5]))
  #check duplicate column
  lst2<- lapply(lst1, function(x) transform(x[!duplicated(x),]))

  ##function to rearrange each df in list in order to count the total coverage
  #in each position and the count of each read nucletide

  all.col= c("seqnames","pos","nucleotide","count"  )
  for (i in seq_along(lst2)){
    colnames(lst2[[i]]) <- all.col
  }
  new= "pos"
  riarrange.df = function(list_df){
    list_df %>%
      dplyr::mutate(end= pos) %>%
      dplyr::group_by(seqnames,pos,end,nucleotide) %>%
      dplyr::summarise(count=sum(count)) %>%
      dplyr::arrange(nucleotide) %>%
      dplyr::summarise(value= paste(sum(count)),
                       counts=paste(nucleotide, count, sep=':',collapse=';'))
  }


  lst3<- lapply(lst2, riarrange.df)

  #pp=Reduce(function(...) merge(...,by= c("seqnames", "pos", "end"), all=TRUE), lst3)
  pp = Reduce(function(x,y) merge(x,y,by=c("seqnames", "pos", "end"),all=TRUE) ,lst3)


  pp[is.na(pp)] <- 0
  pp=as.data.frame(pp)
  if(type_bam == "number"){
    for (i in pp[1]){
      Chromosome<-paste("chr", i, sep="")
      pp<- cbind(Chromosome, pp)
      pp[,2]<-NULL
    }
  }

  ###write file

  dir_users=sprintf("%s/output",outDir)
  myDir <- dir_users
  if (file.exists(myDir)) unlink(myDir,recursive=TRUE)
  dir.create(myDir)
  message("Write output file in directory ", myDir, format(Sys.time(), "%a_%b_%d_%Y"),'.bed')
  message(Sys.time())

  #utils::write.table(x=pp, file =paste0(myDir,format(Sys.time(), "%a_%b_%d_%Y"),'.bed'),
   #           quote=FALSE, sep="\t", eol = "\n", row.names = FALSE,
    #          col.names = FALSE)
  utils::write.table(x=pp, file =file.path(myDir, paste0(format(Sys.time(), "%a_%b_%d_%Y"),'.bed')),
                     quote=FALSE, sep="\t", eol = "\n", row.names = FALSE,
                     col.names = FALSE)

}

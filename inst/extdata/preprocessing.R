#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
output= args[3]
setwd(output)
dir.create('teMpFoldeR', showWarnings= FALSE)

mylist=scan(args[1], character(), quote = "")
#mylist=scan("./mylist.txt", character(), quote = "")
require(AnnotationDbi)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Hs.eg.db)
my_gene_name=select(org.Hs.eg.db, key= mylist, columns="ENTREZID",
                    keytype="ALIAS")


for (i in my_gene_name$ENTREZID){
  if (is.na(i)){
    s=print(subset(my_gene_name$ALIAS,is.na(my_gene_name$ENTREZID)))
    write.table(s, file='./preprocessing_log.txt', quote= FALSE,
                row.names = FALSE, col.names = FALSE)
    stop('ERROR: unrecognized GENE NAME. Please choose a
         correct HGNC official gene names.')}
  else {return= NULL}
}

ID=my_gene_name$ENTREZID
genome= args[2]

if (genome == "hg19"){
  all_gene= TxDb.Hsapiens.UCSC.hg19.knownGene}else{
    all_gene= TxDb.Hsapiens.UCSC.hg38.knownGene}

b=c()
for (i in ID){
  if( ! is.element(i, keys(TxDb.Hsapiens.UCSC.hg19.knownGene, "GENEID")))
    b[i]= i
}
no_entrID= subset(my_gene_name, my_gene_name$ENTREZID %in% b)

if (nrow(no_entrID)!=0){
  write.table(no_entrID, file='./preprocessing_log1.txt',
              quote= FALSE, row.names = FALSE, col.names = FALSE)
  stop('ERROR: unrecognized GENE NAME. Please remove genes
       stored in preprocessing_log1.txt')}


pre= data.frame()
for (i in ID){
  txid <- select(all_gene, i, "TXNAME", "GENEID")[["TXNAME"]]
  cds <-exonsBy(all_gene, by="tx", use.names=TRUE)
  coordinate= as.data.frame(cds[names(cds) %in% txid])
  pre= rbind(coordinate,pre)
  pre$start= pre$start -10
  pre$end= pre$end +10
}


for_bed=data.frame()
for (row in 1:nrow(pre)){
  sel_cc= data.frame(as.character(pre$seqnames[row]),
                     pre$start[row], pre$end[row], stringsAsFactors = FALSE)
  for_bed= rbind(for_bed, sel_cc)
}
colnames(for_bed)[1]='chrm'

type_bam= args[4]

if(type_bam == "number"){
for_bed$chrm= gsub("^.{0,3}", "", for_bed$chrm, perl =TRUE)}

setwd('./teMpFoldeR/')
write.table(x=for_bed, file =paste0(Sys.Date(),'.bed'),
            quote=FALSE, sep="\t", eol = "\n", row.names = FALSE,
            col.names = FALSE)
